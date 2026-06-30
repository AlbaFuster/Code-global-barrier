# MOBIE

MOBIE stands for **Modifying Ocean Barriers for INLA on Earth**.

It is a Shiny app for building and editing global barrier polygons, creating spherical meshes, checking barrier triangles, and exporting outputs for INLA barrier-SPDE workflows.

The app is kept close to the original script. The main code is in:

```text
inst/app/app.R
```

This makes it easier to compare with the working script and avoids changing the behaviour while moving the project to GitHub.

## Run the app

Open the project in RStudio, then run:

```r
devtools::load_all()
run_mobie()
```

You can also run the app directly:

```r
shiny::runApp("inst/app")
```

## Install packages

The app needs these R packages:

```r
install.packages(c(
  "shiny",
  "leaflet",
  "leaflet.extras",
  "sf",
  "rnaturalearth",
  "spdep",
  "units",
  "s2",
  "fmesher",
  "inlabru"
))
```

INLA is installed from the INLA repository:

```r
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

If `INLAspacetime` is not available from your normal repositories, install it from its source repository or the same setup you used for the original script.

## What the app does

- Load a barrier polygon from GeoJSON, JSON, GPKG, or zipped shapefile.
- Build a default world land barrier.
- Simplify and edit the barrier on a Leaflet map.
- Build a spherical mesh.
- Classify barrier triangles.
- Plot local barrier-SPDE correlation diagnostics.
- Export the edited barrier and mesh outputs.

## Project structure

```text
mobie/
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── R/
│   └── run_mobie.R
└── inst/
    └── app/
        └── app.R
```

## Notes

The app uses `fmesher::fmesher_globe_points()` through `library(fmesher)`. This is loaded explicitly because the original app calls the function directly.
