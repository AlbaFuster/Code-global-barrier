# MOBIE

**MOBIE** (**Mo**difying Ocean **B**arriers for **I**NLA on **E**arth) is an R package that provides an interactive Shiny application for creating, editing and exporting global land barriers for barrier SPDE models.

The application was developed to simplify the creation of global barrier geometries and spherical meshes used in spatial analyses with INLA. MOBIE combines interactive editing tools with mesh generation, barrier classification and local correlation diagnostics in a single workflow.

## Main features

MOBIE allows users to:

- Build a default global land barrier from Natural Earth country polygons.
- Import existing barriers from GeoJSON, GeoPackage or zipped shapefiles.
- Edit barrier polygons interactively on a Leaflet map.
- Generate spherical meshes for barrier SPDE models.
- Classify mesh triangles as barrier or non-barrier.
- Explore local correlation patterns of the barrier SPDE model.
- Export edited barriers, meshes and application settings.

## Installation

Install the development version directly from GitHub:

```r
remotes::install_github("AlbaFuster/Code-global-barrier")
```

## Run the application

After installation, launch MOBIE with:

```r
library(mobie)
run_mobie()
```

## Dependencies

MOBIE relies on several R packages for interactive applications, spatial data handling and barrier-SPDE modelling. The main dependencies are:

- **shiny**
- **leaflet** and **leaflet.extras**
- **sf**
- **rnaturalearth**
- **INLA**
- **INLAspacetime**
- **fmesher**
- **inlabru**
- **spdep**
- **s2**
- **units**

## Project structure

```
.
├── R/                  Package functions
├── inst/
│   └── app/            Shiny application
├── DESCRIPTION
├── LICENSE
├── README.md
└── NAMESPACE
```

The Shiny application is located in `inst/app`, while the package provides a simple launcher function (`run_mobie()`) that starts the application after installation.

## Future developments

MOBIE is under active development and new functionality will be added in future releases. Planned improvements include additional tools for barrier editing, mesh construction and support for new spatial modelling workflows.

## Author

**Alba Fuster**

Institut de Ciències del Mar (ICM-CSIC)

## License

This project is distributed under the MIT License.