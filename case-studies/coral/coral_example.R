############################################################
## Global LGCP comparison:
## Barrier SPDE vs standard SPDE for Scleractinia GBIF data
## Version using MOBIE-generated barrier and mesh
############################################################

rm(list = ls())
set.seed(123)

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(fmesher)
library(units)
library(readr)
library(ggpubr)

sf::sf_use_s2(FALSE)

## ---------------------------------------------------------
## 1. CRS definitions
## ---------------------------------------------------------

crs_ll <- st_crs(4326)
crs_vis <- st_crs("+proj=moll +units=km")

## ---------------------------------------------------------
## 2. Load GBIF cleaned data and thin globally
## ---------------------------------------------------------

load("dat_clean_marine.RData")

coral_dat <- dat_clean_marine

coral_dat <- coral_dat %>%
  filter(
    !is.na(decimal_longitude),
    !is.na(decimal_latitude),
    decimal_longitude >= -180,
    decimal_longitude <= 180,
    decimal_latitude >= -90,
    decimal_latitude <= 90
  ) %>%
  distinct(decimal_longitude, decimal_latitude, .keep_all = TRUE)

set.seed(12422)

n_keep <- min(20000, nrow(coral_dat))

coral_dat <- coral_dat %>%
  slice_sample(n = n_keep)

coral_ll <- st_as_sf(
  coral_dat,
  coords = c("decimal_longitude", "decimal_latitude"),
  crs = crs_ll,
  remove = FALSE
)

## ---------------------------------------------------------
## 3. Quick plot of thinned GBIF data
## ---------------------------------------------------------

world_plot <- ne_countries(scale = "medium", returnclass = "sf")

world_mll_plot <- world_plot %>%
  st_transform(crs_vis) %>%
  st_union() %>%
  st_as_sf()

coral_mll <- st_transform(coral_ll, crs_vis)

p_coral <- ggplot() +
  geom_sf(
    data = world_mll_plot,
    fill = "grey85",
    color = "grey45",
    linewidth = 0.15
  ) +
  geom_sf(
    data = coral_mll,
    color = "red",
    size = 0.35,
    alpha = 0.45
  ) +
  coord_sf(crs = crs_vis, expand = FALSE) +
  theme_void() +
  ggtitle("Thinned global Scleractinia GBIF records")

p_coral

## ---------------------------------------------------------
## 4. Load barrier generated with MOBIE
## ---------------------------------------------------------

# In MOBIE, use:
#   Barrier base spherical resolution = 35
#   Simplification tolerance, km = 50
#   Remove small isolated islands = TRUE
#
# Then download:
#   Download barrier GeoJSON

barrier_file <- "barrier_polygon_2026-06-04.geojson"

world_barrier_ll <- st_read(
  barrier_file,
  quiet = TRUE
)

world_barrier_ll <- st_make_valid(world_barrier_ll)

world_barrier <- world_barrier_ll %>%
  st_transform(crs_vis) %>%
  st_geometry() %>%
  st_union() %>%
  st_make_valid()

world_lnd <- world_barrier

## ---------------------------------------------------------
## 5. Load spherical mesh generated with MOBIE
## ---------------------------------------------------------

# In MOBIE, use:
#   Mesh base spherical resolution = 35
#   Mesh globe resolution inside barriers = 10
#   Mesh cutoff multiplier = 1
#
# Then download:
#   Download mesh .RData
#
# The downloaded file should contain:
#   smesh       : spherical mesh object
#   tri_barrier : indices of barrier triangles

load("MOBIE_mesh_2026-06-04.RData")

stopifnot(exists("smesh"))
stopifnot(exists("tri_barrier"))

## ---------------------------------------------------------
## 7. Define barrier and non-barrier SPDE models
## ---------------------------------------------------------

bmodel <- barrierModel.define(
  mesh = smesh,
  barrier.triangles = tri_barrier,
  prior.range = c(0.2, 0.5),
  prior.sigma = c(0.5, 0.5),
  range.fraction = 0.2
)

spde_model <- inla.spde2.pcmatern(
  mesh = smesh,
  alpha = 2,
  prior.range = c(0.2, 0.5),
  prior.sigma = c(0.5, 0.5)
)

## ---------------------------------------------------------
## 8. Prepare point pattern
## ---------------------------------------------------------

coral_sph <- fm_transform(
  coral_ll,
  crs = fm_crs("sphere")
)

coral_sph <- st_sf(
  geometry = st_geometry(coral_sph)
)

## ---------------------------------------------------------
## 9. Global ocean sampler
## ---------------------------------------------------------

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

earth_ll <- Earth_poly(resol = 200)
earth_mll <- st_transform(earth_ll, crs_vis)

ocean_mll <- st_difference(
  earth_mll,
  st_union(world_barrier)
)

ocean_sampler <- st_as_sf(
  st_transform(ocean_mll, crs_ll)
)

ocean_sampler <- st_sf(
  geometry = st_geometry(ocean_sampler)
)

## ---------------------------------------------------------
## 10. Fit LGCP with barrier
## ---------------------------------------------------------

cmp_barrier <- ~
  spatial_barrier(geometry, model = bmodel) +
  Intercept(1)

fit_barrier <- bru(
  cmp_barrier,
  bru_obs(
    formula = geometry ~ .,
    family = "cp",
    data = coral_sph,
    samplers = ocean_sampler,
    domain = list(geometry = smesh)
  ),
  options = list(
    control.inla = list(
      int.strategy = "ccd"
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      return.marginals.predictor = TRUE
    )
  )
)

summary(fit_barrier)

## ---------------------------------------------------------
## 11. Fit LGCP without barrier
## ---------------------------------------------------------

cmp_nobarrier <- ~
  spatial_nobarrier(geometry, model = spde_model) +
  Intercept(1)

fit_nobarrier <- bru(
  cmp_nobarrier,
  bru_obs(
    formula = geometry ~ .,
    family = "cp",
    data = coral_sph,
    samplers = ocean_sampler,
    domain = list(geometry = smesh)
  ),
  options = list(
    control.inla = list(
      int.strategy = "ccd"
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      return.marginals.predictor = TRUE
    )
  )
)

summary(fit_nobarrier)

## ---------------------------------------------------------
## 12. Model comparison table
## ---------------------------------------------------------

model_comparison <- data.frame(
  model = c("Barrier SPDE", "Standard SPDE"),
  waic = c(
    fit_barrier$waic$waic,
    fit_nobarrier$waic$waic
  ),
  dic = c(
    fit_barrier$dic$dic,
    fit_nobarrier$dic$dic
  ),
  marginal_loglik = c(
    fit_barrier$mlik[1],
    fit_nobarrier$mlik[1]
  )
)

print(model_comparison)

## ---------------------------------------------------------
## 13. Prediction grid in Mollweide
## ---------------------------------------------------------

grid_create <- function(barrier, resol = 150) {
  
  Ell <- Earth_poly(resol = 100)
  Ell_mll <- st_transform(Ell, st_crs(barrier))
  
  xy <- expand.grid(
    x = seq(
      st_bbox(Ell_mll)["xmin"],
      st_bbox(Ell_mll)["xmax"],
      by = resol
    ),
    y = seq(
      st_bbox(Ell_mll)["ymin"],
      st_bbox(Ell_mll)["ymax"],
      by = resol
    )
  )
  
  grid <- st_as_sf(
    xy,
    coords = c("x", "y"),
    crs = st_crs(barrier)
  )
  
  inside_earth <- lengths(st_intersects(grid, Ell_mll)) > 0
  inside_barrier <- lengths(st_intersects(grid, barrier)) > 0
  
  grid[inside_earth & !inside_barrier, ]
}

pred_grid_mll <- grid_create(
  barrier = world_barrier,
  resol = 150
)

pred_grid_ll <- st_transform(pred_grid_mll, crs_ll)

pred_grid_sph <- fm_transform(
  pred_grid_ll,
  crs = fm_crs("sphere")
)

pred_grid_sph <- st_sf(
  geometry = st_geometry(pred_grid_sph)
)

## ---------------------------------------------------------
## 14. Predictions: log intensity and spatial effects
## ---------------------------------------------------------

pred_barrier <- predict(
  fit_barrier,
  newdata = pred_grid_sph,
  formula = ~ Intercept + spatial_barrier
)

pred_nobarrier <- predict(
  fit_nobarrier,
  newdata = pred_grid_sph,
  formula = ~ Intercept + spatial_nobarrier
)

pred_spatial_barrier <- predict(
  fit_barrier,
  newdata = pred_grid_sph,
  formula = ~ spatial_barrier
)

pred_spatial_nobarrier <- predict(
  fit_nobarrier,
  newdata = pred_grid_sph,
  formula = ~ spatial_nobarrier
)

pred_grid_mll$logint_barrier <- pred_barrier$mean
pred_grid_mll$logint_barrier_sd <- pred_barrier$sd
pred_grid_mll$logint_nobarrier <- pred_nobarrier$mean
pred_grid_mll$logint_nobarrier_sd <- pred_nobarrier$sd

pred_grid_mll$spatial_barrier <- pred_spatial_barrier$mean
pred_grid_mll$spatial_nobarrier <- pred_spatial_nobarrier$mean

pred_grid_mll$diff_logint <- pred_grid_mll$logint_barrier -
  pred_grid_mll$logint_nobarrier

pred_grid_mll$diff_spatial <- pred_grid_mll$spatial_barrier -
  pred_grid_mll$spatial_nobarrier

## Scale log-intensity maps jointly to 0-1 for visual comparison

joint_rng <- range(
  c(
    pred_grid_mll$logint_barrier,
    pred_grid_mll$logint_nobarrier
  ),
  na.rm = TRUE
)

pred_grid_mll$logint_barrier_01 <-
  (pred_grid_mll$logint_barrier - joint_rng[1]) /
  (joint_rng[2] - joint_rng[1])

pred_grid_mll$logint_nobarrier_01 <-
  (pred_grid_mll$logint_nobarrier - joint_rng[1]) /
  (joint_rng[2] - joint_rng[1])

pred_grid_mll$diff_logint_01 <- pred_grid_mll$logint_barrier_01 -
  pred_grid_mll$logint_nobarrier_01

## ---------------------------------------------------------
## 15. Prepare plotting data
## ---------------------------------------------------------

pred_df <- pred_grid_mll %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

tile_width <- 150
tile_height <- 150

## ---------------------------------------------------------
## 16. Plots: barrier vs no barrier
## ---------------------------------------------------------

p_barrier <- ggplot() +
  geom_tile(
    data = pred_df,
    aes(x = x, y = y, fill = logint_barrier_01),
    width = tile_width,
    height = tile_height
  ) +
  geom_sf(
    data = world_mll_plot,
    fill = "grey85",
    color = "grey45",
    linewidth = 0.15
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    limits = c(0, 1),
    name = "Scaled\nlog-intensity"
  ) +
  coord_sf(crs = crs_vis, expand = FALSE) +
  theme_void() +
  ggtitle("B) Barrier model")

p_nobarrier <- ggplot() +
  geom_tile(
    data = pred_df,
    aes(x = x, y = y, fill = logint_nobarrier_01),
    width = tile_width,
    height = tile_height
  ) +
  geom_sf(
    data = world_mll_plot,
    fill = "grey85",
    color = "grey45",
    linewidth = 0.15
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    limits = c(0, 1),
    name = "Scaled\nlog-intensity"
  ) +
  coord_sf(crs = crs_vis, expand = FALSE) +
  theme_void() +
  ggtitle("C) Standard SPDE model")

diff_lim <- max(abs(pred_df$diff_logint_01), na.rm = TRUE)

p_diff <- ggplot() +
  geom_tile(
    data = pred_df,
    aes(x = x, y = y, fill = diff_logint_01),
    width = tile_width,
    height = tile_height
  ) +
  geom_sf(
    data = world_mll_plot,
    fill = "grey85",
    color = "grey45",
    linewidth = 0.15
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-diff_lim, diff_lim),
    name = "Barrier -\nstandard"
  ) +
  coord_sf(crs = crs_vis, expand = FALSE) +
  theme_void() +
  ggtitle("D) Difference in scaled log-intensity")

p_compare <- ggarrange(
  p_coral,
  p_barrier,
  p_nobarrier,
  p_diff,
  nrow = 2,
  ncol = 2
)

p_compare

## ---------------------------------------------------------
## 17. Plots: spatial effects
## ---------------------------------------------------------

sp_lim <- max(
  abs(c(pred_df$spatial_barrier, pred_df$spatial_nobarrier)),
  na.rm = TRUE
)

p_spatial_barrier <- ggplot() +
  geom_tile(
    data = pred_df,
    aes(x = x, y = y, fill = spatial_barrier),
    width = tile_width,
    height = tile_height
  ) +
  geom_sf(
    data = world_mll_plot,
    fill = "grey85",
    color = "grey45",
    linewidth = 0.15
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-sp_lim, sp_lim),
    name = "Spatial\neffect"
  ) +
  coord_sf(crs = crs_vis, expand = FALSE) +
  theme_void() +
  ggtitle("Barrier spatial effect")

p_spatial_nobarrier <- ggplot() +
  geom_tile(
    data = pred_df,
    aes(x = x, y = y, fill = spatial_nobarrier),
    width = tile_width,
    height = tile_height
  ) +
  geom_sf(
    data = world_mll_plot,
    fill = "grey85",
    color = "grey45",
    linewidth = 0.15
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-sp_lim, sp_lim),
    name = "Spatial\neffect"
  ) +
  coord_sf(crs = crs_vis, expand = FALSE) +
  theme_void() +
  ggtitle("Standard SPDE spatial effect")

p_spatial_diff <- ggplot() +
  geom_tile(
    data = pred_df,
    aes(x = x, y = y, fill = diff_spatial),
    width = tile_width,
    height = tile_height
  ) +
  geom_sf(
    data = world_mll_plot,
    fill = "grey85",
    color = "grey45",
    linewidth = 0.15
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = "Barrier -\nstandard"
  ) +
  coord_sf(crs = crs_vis, expand = FALSE) +
  theme_void() +
  ggtitle("Difference in spatial effect")

p_spatial_compare <- ggarrange(
  p_spatial_barrier,
  p_spatial_nobarrier,
  p_spatial_diff,
  nrow = 1,
  ncol = 3
)

p_spatial_compare