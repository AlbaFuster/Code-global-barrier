############################################################
## HURDLE SPATIAL MODEL FOR Gadus morhua
## Presence/absence + positive biomass
## Version using MOBIE-generated barrier and mesh
############################################################

rm(list = ls())

set.seed(52134)

## ============================================================
## 0. LIBRARIES
## ============================================================

library(sf)
library(sp)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(INLA)
library(INLAspacetime)
library(fmesher)
library(spdep)
library(viridis)
library(patchwork)

sf::sf_use_s2(FALSE)

## ============================================================
## 1. CRS DEFINITIONS AND MOBIE FILES
## ============================================================

# Longitude-latitude CRS.
crs_ll <- st_crs(4326)

# Mollweide projection in kilometres.
# This is used for global plotting, distance operations, and grid construction.
crs_vis <- st_crs("+proj=moll +units=km")

# Files exported from MOBIE.
#
# MOBIE parameters used to reproduce the original script:
#
# Barrier:
#   Barrier base spherical resolution = 50
#   Simplification tolerance, km      = 20
#   Remove small isolated islands     = TRUE
#
# Mesh:
#   Mesh base spherical resolution       = 50
#   Mesh globe resolution inside barriers = 30
#   Mesh cutoff multiplier                = 1

barrier_file <- "barrier_polygon_2026-06-04.geojson"
mesh_file <- "MOBIE_mesh_2026-06-04.RData"

## ============================================================
## 2. LOAD AND PREPARE FISHGLOB DATA
## ============================================================

# Load FishGlob public data.
data <- read.csv("FishGlob_public_clean.csv")

# Keep the relevant variables and enforce consistent types.
data_specific <- data.frame(
  species   = as.character(data$accepted_name),
  longitude = as.numeric(data$longitude),
  latitude  = as.numeric(data$latitude),
  CPUA      = as.numeric(data$num_cpua),
  CPUE      = as.numeric(data$num_cpue),
  year      = as.numeric(data$year),
  month     = as.numeric(data$month),
  survey    = as.character(data$survey),
  country   = as.character(data$country)
)

# Remove rows with missing values.
data_specific <- na.omit(data_specific)

# Select Atlantic cod records from 2000 onwards.
data_species <- data_specific %>%
  filter(
    species == "Gadus morhua",
    year >= 2000,
    !is.na(CPUA),
    !is.na(longitude),
    !is.na(latitude)
  )

# Convert observations to an sf object in longitude-latitude coordinates.
# Duplicate spatial locations are removed to avoid repeated identical geometries.
data_species_sf <- st_as_sf(
  data_species,
  coords = c("longitude", "latitude"),
  crs = crs_ll
) %>%
  distinct(geometry, .keep_all = TRUE)

# Positive biomass observations for the lognormal part of the hurdle model.
data_pos_sf <- data_species_sf %>%
  filter(CPUA > 0)

# Presence locations for the Bernoulli presence/absence model.
presence_sf <- data_pos_sf %>%
  mutate(pa = 1)

## ============================================================
## 3. LOAD MOBIE BARRIER AND CREATE OCEAN GRID
## ============================================================

# The MOBIE barrier is exported in longitude-latitude coordinates.
world_barrier_ll <- st_read(
  barrier_file,
  quiet = TRUE
)

# Validate the geometry to avoid topology issues.
world_barrier_ll <- st_make_valid(world_barrier_ll)

# Transform the barrier to Mollweide because the prediction grid and
# absence selection are created in projected space.
world_barrier <- world_barrier_ll %>%
  st_transform(crs_vis) %>%
  st_geometry() %>%
  st_union() %>%
  st_make_valid()

# Keep this object name for consistency with the original script.
world_lnd <- world_barrier

# Land layer used in the plots.
world_mll <- world_barrier

# Function defining an approximate global polygon.
# This is used to create an ocean domain by subtracting the land barrier.
Earth_poly <- function(resol = 100) {
  st_sfc(
    st_multipolygon(
      list(
        st_polygon(
          list(
            cbind(
              long = c(
                seq(-1, 1, length.out = resol * 2 + 1),
                rep(1, resol + 1),
                seq(1, -1, length.out = resol * 2 + 1),
                rep(-1, resol + 1)
              ) * (180 - 1e-5),
              lat = c(
                rep(1, resol * 2 + 1),
                seq(1, -1, length.out = resol + 1),
                rep(-1, resol * 2 + 1),
                seq(-1, 1, length.out = resol + 1)
              ) * (90 - 1e-5)
            )
          )
        )
      )
    ),
    crs = crs_ll
  )
}

# Create a regular prediction grid over the ocean.
# The grid is first created in Mollweide coordinates and then land cells are removed.
grid_create <- function(barrier, resol = 50) {
  
  Ell <- Earth_poly(resol = 100)
  Ell <- st_transform(Ell, st_crs(barrier))
  
  grid_0 <- st_as_sf(
    expand.grid(
      x = seq(
        min(st_coordinates(Ell)[, 1]),
        max(st_coordinates(Ell)[, 1]),
        by = resol
      ),
      y = seq(
        min(st_coordinates(Ell)[, 2]),
        max(st_coordinates(Ell)[, 2]),
        by = resol
      )
    ),
    coords = c("x", "y"),
    crs = st_crs(barrier)
  )
  
  inside_earth <- lengths(st_intersects(grid_0, Ell)) > 0
  on_land <- lengths(st_intersects(grid_0, barrier)) > 0
  
  grid_0[inside_earth & !on_land, ]
}

# Prediction grid in Mollweide.
mgrid <- grid_create(
  barrier = world_barrier,
  resol = 50
)

# Same grid in longitude-latitude coordinates.
grid_ll <- st_transform(
  mgrid,
  crs_ll
)

# Same grid transformed to the sphere for INLA/fmesher projection.
mgrid_sph <- fm_transform(
  mgrid,
  crs = fm_crs("sphere")
)

## ============================================================
## 4. LOAD MOBIE SPHERICAL MESH AND DEFINE BARRIER MODEL
## ============================================================

# MOBIE exports an .RData file containing:
#   smesh       : spherical mesh object
#   tri_barrier : vector of barrier triangle indices

load(mesh_file)

stopifnot(exists("smesh"))
stopifnot(exists("tri_barrier"))

# Keep the original object name used in the bacalao script.
triBarrier <- tri_barrier

# Define the barrier SPDE model.
# The mesh and barrier triangles come directly from MOBIE.
bmodel <- barrierModel.define(
  mesh = smesh,
  barrier.triangles = triBarrier,
  prior.range = c(0.5, 0.5),
  prior.sigma = c(1, 0.5),
  range.fraction = 0.008,
  constr = TRUE
)

# Projector from mesh nodes to prediction grid.
gproj <- inla.mesh.projector(
  mesh = smesh,
  loc = mgrid_sph
)

## ============================================================
## 5. POSITIVE BIOMASS MODEL: LOGNORMAL COMPONENT
## ============================================================

# Positive-biomass locations transformed to the sphere.
locs_pos_sph <- fm_transform(
  data_pos_sf$geometry,
  crs = fm_crs("sphere")
)

# Projection matrix from mesh nodes to positive-biomass observations.
A_pos <- inla.spde.make.A(
  mesh = smesh,
  loc = locs_pos_sph
)

data_lognormal <- data.frame(
  y = data_pos_sf$CPUA
)

# Estimation stack for positive biomass.
stk_pos <- inla.stack(
  data = list(y = data_lognormal$y),
  A = list(A_pos, 1),
  effects = list(
    i_biomass = 1:bmodel$mesh$n,
    beta0_biomass = rep(1, nrow(data_lognormal))
  ),
  tag = "est_biomass"
)

# Prediction stack for positive biomass over the ocean grid.
stk_pred_biomass <- inla.stack(
  data = list(y = NA),
  A = list(gproj$proj$A, 1),
  effects = list(
    i_biomass = 1:bmodel$mesh$n,
    beta0_biomass = rep(1, nrow(mgrid))
  ),
  tag = "pred_biomass"
)

stk_full_biomass <- inla.stack(
  stk_pos,
  stk_pred_biomass
)

# Fit lognormal model to positive CPUA values.
res_biomass <- inla(
  y ~ -1 + beta0_biomass +
    f(i_biomass, model = bmodel),
  data = inla.stack.data(stk_full_biomass),
  family = "lognormal",
  control.predictor = list(
    A = inla.stack.A(stk_full_biomass),
    link = 1,
    compute = TRUE
  ),
  control.compute = list(
    config = TRUE,
    return.marginals.predictor = TRUE
  ),
  verbose = TRUE
)

## ============================================================
## 6. GENERATE ABSENCES FAR FROM PRESENCES
## ============================================================

# Transform presences and grid to Mollweide for distance calculations.
presence_mll <- st_transform(
  presence_sf,
  crs_vis
)

grid_mll <- st_transform(
  grid_ll,
  crs_vis
)

# Minimum distance from presences, in kilometres.
min_dist_presence <- 2000

# Candidate absence grid cells must be far from observed presences.
far_from_presence <- lengths(
  st_is_within_distance(
    grid_mll,
    presence_mll,
    dist = min_dist_presence
  )
) == 0

candidate_absences <- grid_mll[far_from_presence, ]

# Use up to three times as many absences as presences.
n_pres <- nrow(presence_sf)
n_abs <- min(
  nrow(candidate_absences),
  n_pres * 3
)

pseudo_abs_mll <- candidate_absences[
  sample(seq_len(nrow(candidate_absences)), n_abs),
]

pseudo_abs_sf <- st_transform(
  pseudo_abs_mll,
  crs_ll
) %>%
  mutate(pa = 0)

presence_pa_sf <- presence_sf %>%
  select(pa, geometry)

pa_sf <- rbind(
  presence_pa_sf,
  pseudo_abs_sf[, c("pa", "geometry")]
)

## ============================================================
## 7. PRESENCE/ABSENCE MODEL: BINOMIAL COMPONENT
## ============================================================

# Presence/absence locations transformed to the sphere.
locs_pa_sph <- fm_transform(
  pa_sf$geometry,
  crs = fm_crs("sphere")
)

# Projection matrix from mesh nodes to presence/absence observations.
A_pa <- inla.spde.make.A(
  mesh = smesh,
  loc = locs_pa_sph
)

data_pa <- data.frame(
  y = pa_sf$pa
)

# Estimation stack for Bernoulli model.
stk_pa <- inla.stack(
  data = list(y = data_pa$y),
  A = list(A_pa, 1),
  effects = list(
    i_pa = 1:bmodel$mesh$n,
    beta0_pa = rep(1, nrow(data_pa))
  ),
  tag = "est_pa"
)

# Prediction stack for presence probability over the ocean grid.
stk_pred_pa <- inla.stack(
  data = list(y = NA),
  A = list(gproj$proj$A, 1),
  effects = list(
    i_pa = 1:bmodel$mesh$n,
    beta0_pa = rep(1, nrow(mgrid))
  ),
  tag = "pred_pa"
)

stk_full_pa <- inla.stack(
  stk_pa,
  stk_pred_pa
)

# Fit Bernoulli model to presences and absences.
res_pa <- inla(
  y ~ -1 + beta0_pa +
    f(i_pa, model = bmodel),
  data = inla.stack.data(stk_full_pa),
  family = "binomial",
  control.predictor = list(
    A = inla.stack.A(stk_full_pa),
    link = 1,
    compute = TRUE
  ),
  control.compute = list(
    config = TRUE,
    return.marginals.predictor = TRUE
  ),
  verbose = TRUE
)

## ============================================================
## 8. EXTRACT PREDICTIONS AND HURDLE CORRECTION
## ============================================================

# Prediction indices in each INLA stack.
idx_biomass <- inla.stack.index(
  stk_full_biomass,
  "pred_biomass"
)$data

idx_pa <- inla.stack.index(
  stk_full_pa,
  "pred_pa"
)$data

# Positive biomass predictions.
pred_biomass <- res_biomass$summary.fitted.values$mean[idx_biomass]
pred_biomass_sd <- res_biomass$summary.fitted.values$sd[idx_biomass]

# Presence probability predictions.
pred_prob <- res_pa$summary.fitted.values$mean[idx_pa]
pred_prob_sd <- res_pa$summary.fitted.values$sd[idx_pa]

# Hurdle-corrected biomass:
# expected biomass = positive biomass conditional on presence * probability of presence.
pred_hurdle <- pred_biomass * pred_prob

data_pred <- data.frame(
  x = st_coordinates(grid_ll)[, 1],
  y = st_coordinates(grid_ll)[, 2],
  biomass_positive = pred_biomass,
  biomass_positive_sd = pred_biomass_sd,
  prob_presence = pred_prob,
  prob_presence_sd = pred_prob_sd,
  biomass_hurdle = pred_hurdle
)

data_pred_sf <- st_as_sf(
  data_pred,
  coords = c("x", "y"),
  crs = crs_ll
)

data_pred_mll <- st_transform(
  data_pred_sf,
  crs_vis
)

## ============================================================
## 9. PLOTS
## ============================================================

# Soft blue-grey palette for biomass.
pal_biomass <- c(
  "#d9e3ea",
  "#a6bddb",
  "#6b8fb3",
  "#3e5f7d"
)

# Neutral green palette for occurrence probability.
pal_prob <- c(
  "#e5eadf",
  "#c6d4b8",
  "#91ad7e",
  "#5f7f55"
)

land_layer <- geom_sf(
  data = world_mll,
  fill = "grey85",
  color = "grey60",
  alpha = 0.9
)

theme_map <- theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 11
    ),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# A) Predicted relative biomass conditional on presence.
p_a <- ggplot() +
  theme_map +
  geom_sf(
    aes(color = biomass_positive),
    data = data_pred_mll,
    size = 0.7
  ) +
  scale_color_gradientn(
    colours = pal_biomass,
    name = "Relative\nbiomass"
  ) +
  land_layer +
  ggtitle("A) Predicted relative biomass")

# B) Predicted occurrence probability.
p_b <- ggplot() +
  theme_map +
  geom_sf(
    aes(color = prob_presence),
    data = data_pred_mll,
    size = 0.7
  ) +
  scale_color_gradientn(
    colours = pal_prob,
    limits = c(0, 1),
    name = "Presence\nprobability"
  ) +
  land_layer +
  ggtitle("B) Predicted probability")

# C) Hurdle-corrected relative biomass.
p_c <- ggplot() +
  theme_map +
  geom_sf(
    aes(color = biomass_hurdle),
    data = data_pred_mll,
    size = 0.7
  ) +
  scale_color_gradientn(
    colours = pal_biomass,
    name = "Corrected\nbiomass"
  ) +
  land_layer +
  ggtitle("C) Hurdle-corrected relative biomass")

# Final panel:
# first row contains biomass and probability;
# second row contains hurdle-corrected biomass across both columns.
panel_hurdle <- (p_a | p_b) / p_c +
  plot_layout(
    heights = c(1, 1.15)
  )

panel_hurdle