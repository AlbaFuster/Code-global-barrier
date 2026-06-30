############################################################
## Download and clean GBIF records for Scleractinia
############################################################

## Packages
pkgs <- c(
  "rgbif", "dplyr", "ggplot2", "sf", "rnaturalearth",
  "CoordinateCleaner", "readr"
)

install.packages(setdiff(pkgs, rownames(installed.packages())))

library(rgbif)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(CoordinateCleaner)
library(readr)

## ---------------------------------------------------------
## 1. GBIF credentials
## ---------------------------------------------------------
## Create an account at https://www.gbif.org/
## Then replace these values:
user <- "YOUR USER"
pwd  <- "YOUR PASSEWORD"
email <- "YOUR EMAIL"

## ---------------------------------------------------------
## 2. Get taxon key for Scleractinia
## ---------------------------------------------------------
tax <- name_backbone(name = "Scleractinia", rank = "ORDER")
taxon_key <- tax$usageKey

taxon_key

## ---------------------------------------------------------
## 3. Download GBIF occurrences
## ---------------------------------------------------------
## Filters applied directly in GBIF:
## - only preserved/published occurrence records
## - only records with coordinates
## - remove records already flagged as geospatially problematic
## - remove records with zero coordinate uncertainty if suspicious
## - remove fossils/living specimens if not desired

download_key <- occ_download(
  pred("taxonKey", taxon_key),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_in("basisOfRecord", c(
    "HUMAN_OBSERVATION",
    "OBSERVATION",
    "MACHINE_OBSERVATION",
    "PRESERVED_SPECIMEN",
    "MATERIAL_SAMPLE"
  )),
  pred("occurrenceStatus", "PRESENT"),
  format = "SIMPLE_CSV",
  user = user,
  pwd = pwd,
  email = email
)

## Wait until GBIF prepares the download
occ_download_wait(download_key)

## Download and import
dir.create("data_gbif", showWarnings = FALSE)

gbif_zip <- occ_download_get(download_key, path = "data_gbif/")

dat_raw <- occ_download_import(gbif_zip, path = "./data_gbif/")

## ---------------------------------------------------------
## 4. Basic cleaning
## ---------------------------------------------------------
## ---------------------------------------------------------
## 4. Basic cleaning + tropical coral regions
## ---------------------------------------------------------

dat <- dat_raw %>%
  janitor::clean_names() %>%
  
  ## Valid coordinates
  filter(
    !is.na(decimal_longitude),
    !is.na(decimal_latitude),
    decimal_longitude >= -180,
    decimal_longitude <= 180,
    decimal_latitude >= -90,
    decimal_latitude <= 90
  ) %>%
  
  ## Coordinate uncertainty
  filter(
    is.na(coordinate_uncertainty_in_meters) |
      coordinate_uncertainty_in_meters <= 50000
  ) %>%
  
  ## Recent observations
  filter(
    is.na(year) | year >= 1950
  ) %>%
  
  ## Tropical/subtropical coral belt
  filter(
    
    ## Indo-Pacific / Coral Triangle
    (
      decimal_longitude >= 90 &
        decimal_longitude <= 180 &
        decimal_latitude >= -35 &
        decimal_latitude <= 35
    ) |
      
      ## Red Sea + Arabian Sea
      (
        decimal_longitude >= 30 &
          decimal_longitude <= 75 &
          decimal_latitude >= -5 &
          decimal_latitude <= 30
      ) |
      
      ## Western Indian Ocean / Madagascar
      (
        decimal_longitude >= 35 &
          decimal_longitude <= 65 &
          decimal_latitude >= -30 &
          decimal_latitude <= 5
      ) |
      
      ## Tropical eastern Pacific
      (
        decimal_longitude >= -120 &
          decimal_longitude <= -75 &
          decimal_latitude >= -10 &
          decimal_latitude <= 30
      ) |
      
      ## Caribbean
      (
        decimal_longitude >= -90 &
          decimal_longitude <= -55 &
          decimal_latitude >= 5 &
          decimal_latitude <= 30
      ) |
      
      ## Brazil tropical reefs
      (
        decimal_longitude >= -50 &
          decimal_longitude <= -30 &
          decimal_latitude >= -25 &
          decimal_latitude <= 5
      )
    
  ) %>%
  
  ## Remove duplicated observations
  distinct(
    species,
    decimal_longitude,
    decimal_latitude,
    year,
    month,
    day,
    .keep_all = TRUE
  )

## ---------------------------------------------------------
## 5. CoordinateCleaner filters
## ---------------------------------------------------------
## Flags suspicious coordinates:
## - zeros
## - country centroids
## - capitals
## - GBIF headquarters
## - biodiversity institutions
## - invalid coordinates
## - equal lat/lon

flags <- clean_coordinates(
  x = dat,
  lon = "decimal_longitude",
  lat = "decimal_latitude",
  species = "species",
  tests = c(
    "capitals",
    "centroids",
    "equal",
    "gbif",
    "institutions",
    "zeros"
  ),
  value = "spatialvalid"
)

dat_clean <- dat[flags$.summary, ]

## ---------------------------------------------------------
## 6. Optional: keep only marine records
## ---------------------------------------------------------
## Scleractinia are marine, so remove records falling on land.
## This is useful because GBIF sometimes contains imprecise coordinates.

world <- ne_countries(scale = "medium", returnclass = "sf")

pts <- st_as_sf(
  dat_clean,
  coords = c("decimal_longitude", "decimal_latitude"),
  crs = 4326,
  remove = FALSE
)

on_land <- lengths(st_intersects(pts, world)) > 0

dat_clean_marine <- dat_clean[!on_land, ]

## ---------------------------------------------------------
## 7. Save outputs
## ---------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)

write_csv(dat_raw, "outputs/scleractinia_gbif_raw.csv")
write_csv(dat_clean, "outputs/scleractinia_gbif_clean.csv")
write_csv(dat_clean_marine, "outputs/scleractinia_gbif_clean_marine.csv")

save(
  dat_raw,
  dat_clean,
  dat_clean_marine,
  file = "outputs/scleractinia_gbif_cleaned.RData"
)

## ---------------------------------------------------------
## 8. Plot
## ---------------------------------------------------------
world_plot <- ne_countries(scale = "medium", returnclass = "sf")

p_raw <- ggplot() +
  geom_sf(data = world_plot, fill = "grey85", color = "grey50", linewidth = 0.2) +
  geom_point(
    data = dat,
    aes(decimal_longitude, decimal_latitude),
    size = 0.4,
    alpha = 0.3
  ) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("Scleractinia GBIF records before coordinate cleaning")

p_clean <- ggplot() +
  geom_sf(data = world_plot, fill = "grey85", color = "grey50", linewidth = 0.2) +
  geom_point(
    data = dat_clean_marine,
    aes(decimal_longitude, decimal_latitude),
    size = 0.4,
    alpha = 0.5
  ) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("Scleractinia GBIF records after cleaning")

p_raw
p_clean

ggsave("outputs/scleractinia_raw_map.png", p_raw, width = 10, height = 5, dpi = 300)
ggsave("outputs/scleractinia_clean_marine_map.png", p_clean, width = 10, height = 5, dpi = 300)

## ---------------------------------------------------------
## 9. Summary
## ---------------------------------------------------------
cat("Raw records: ", nrow(dat_raw), "\n")
cat("After basic filters: ", nrow(dat), "\n")
cat("After CoordinateCleaner: ", nrow(dat_clean), "\n")
cat("After removing land records: ", nrow(dat_clean_marine), "\n")