# ==============================================================================
# MOBIE: Modifying Ocean Barriers for INLA on Earth
# ==============================================================================
#
# Purpose:
#   This Shiny app supports the construction, manual editing, visualization,
#   export, and diagnostic evaluation of global ocean/land barrier polygons for
#   barrier SPDE models on the sphere.
#
# Main functionality:
#   1. Load an external polygon barrier from GeoJSON, JSON, GPKG, or ZIP shapefile.
#   2. Build a default global land barrier from rnaturalearth country polygons.
#   3. Simplify and optionally remove small isolated islands from the barrier.
#   4. Draw polygons interactively on a Leaflet map.
#   5. Add drawn polygons to the barrier or subtract them from it.
#   6. Build a spherical INLA/fmesher mesh over the ocean domain.
#   7. Classify mesh triangles as barrier or non-barrier triangles.
#   8. Compute local barrier-SPDE correlation diagnostics from clicked locations.
#   9. Export the edited barrier as GeoJSON.
#  10. Export the spherical mesh and barrier-triangle index as an RData file.
#
# Notes:
#   - All user-facing documentation is written in English.
#   - Geometry editing is done in Mollweide projection because distances and
#     simplification tolerances are expressed in kilometres.
#   - Spherical mesh construction uses fmesher/INLA-compatible objects.
#   - The app disables sf's s2 geometry engine because several planar operations
#     such as st_difference(), st_union(), st_simplify(), and buffering are more
#     predictable here in a projected CRS.
#
# ==============================================================================


# ==============================================================================
# 1. Libraries
# ==============================================================================

# Core Shiny framework for building the web application.
library(shiny)

# Interactive web maps.
library(leaflet)

# Drawing tools for Leaflet, used to let the user draw polygons.
library(leaflet.extras)

# Simple Features package for reading, transforming, validating, and manipulating
# geospatial vector data.
library(sf)

# Provides Natural Earth country polygons used to create the default global barrier.
library(rnaturalearth)

# Spatial neighbourhood tools; used here to detect isolated small land polygons.
library(spdep)

# File extension utilities used when reading uploaded files.
library(tools)

# Unit handling for spatial quantities such as polygon area.
library(units)

# Spherical geometry library used for nearest-neighbour distance calculations
# on the Earth sphere in the local correlation diagnostic.
library(s2)

# INLA is required for SPDE/barrier-SPDE mesh tools and sparse precision methods.
library(INLA)

# INLAspacetime provides barrier-related tools used by the mesh/SPDE workflow.
library(INLAspacetime)

# fmesher provides mesh utilities used directly by the app, including
# fmesher_globe_points().
library(fmesher)

# inlabru/fmesher provides modern mesh utilities, CRS handling, mesh projection,
# and spherical mesh construction helpers.
library(inlabru)


# ==============================================================================
# 2. Global options and constants
# ==============================================================================

# Allow large uploads, for example zipped shapefiles or high-resolution GeoJSONs.
options(shiny.maxRequestSize = 100 * 1024^2)

# Disable s2 within sf. This app performs several planar operations after
# projection to Mollweide, so GEOS-based planar operations are preferred.
sf::sf_use_s2(FALSE)

# Default spherical resolution used when constructing the default land barrier.
DEFAULT_BARRIER_SRESOL <- 30

# Default spherical resolution used when creating candidate mesh locations.
DEFAULT_MESH_SRESOL <- 30

# Mollweide CRS in kilometres. This projection is useful for global visualization
# and approximate distance-based operations such as simplification and buffering.
CRS_MOLL_KM <- st_crs("+proj=moll +units=km")

# Geographic longitude/latitude CRS.
CRS_LONG_LAT <- 4326

# Approximate Earth equatorial circumference in kilometres.
# Used only to define resolution-dependent approximate mesh/barrier scales.
EARTH_CIRCUMFERENCE_KM <- 36080.2


# ==============================================================================
# 3. General geometry helper functions
# ==============================================================================

# ------------------------------------------------------------------------------
# validate_geometry()
# ------------------------------------------------------------------------------
# Purpose:
#   Ensures that an sf/sfc object has valid geometries.
#
# Why:
#   Uploaded, simplified, unioned, or differenced polygons can contain invalid
#   rings, self-intersections, or topology problems. Calling st_make_valid()
#   reduces the risk of downstream failures.
# ------------------------------------------------------------------------------
validate_geometry <- function(x) {
  st_make_valid(x)
}


# ------------------------------------------------------------------------------
# transform_to_longlat()
# ------------------------------------------------------------------------------
# Purpose:
#   Convert an sf/sfc object to EPSG:4326 longitude/latitude coordinates.
#
# Details:
#   If the object has no CRS, EPSG:4326 is assumed. This is useful for user uploads,
#   but it is also a possible source of error if the uploaded file uses another CRS
#   without declaring it.
# ------------------------------------------------------------------------------
transform_to_longlat <- function(x) {
  if (is.na(st_crs(x))) {
    warning("Input geometry has no CRS. Assuming EPSG:4326.")
    st_crs(x) <- CRS_LONG_LAT
  }
  
  validate_geometry(st_transform(validate_geometry(x), CRS_LONG_LAT))
}


# ------------------------------------------------------------------------------
# transform_to_mollweide()
# ------------------------------------------------------------------------------
# Purpose:
#   Convert an sf/sfc object to the Mollweide CRS in kilometres.
#
# Why:
#   Barrier editing, simplification, and approximate buffering are more stable in
#   a projected CRS than in raw longitude/latitude coordinates.
# ------------------------------------------------------------------------------
transform_to_mollweide <- function(x) {
  validate_geometry(st_transform(validate_geometry(x), CRS_MOLL_KM))
}


# ------------------------------------------------------------------------------
# union_geometries()
# ------------------------------------------------------------------------------
# Purpose:
#   Safely union sf/sfc geometries into a single geometry.
#
# Handles:
#   - A single sf object.
#   - A single sfc object.
#   - A list of sf/sfc objects.
#
# Returns:
#   A valid sfc geometry.
# ------------------------------------------------------------------------------
union_geometries <- function(x) {
  if (is.list(x) && !inherits(x, c("sf", "sfc"))) {
    geoms <- do.call(
      c,
      lapply(x, function(z) st_geometry(validate_geometry(z)))
    )
    
    return(validate_geometry(st_union(geoms)))
  }
  
  validate_geometry(st_union(st_geometry(validate_geometry(x))))
}


# ------------------------------------------------------------------------------
# difference_geometries()
# ------------------------------------------------------------------------------
# Purpose:
#   Subtract geometry y from geometry x.
#
# Why:
#   Used when the user draws a polygon and clicks "Cut drawn polygon".
# ------------------------------------------------------------------------------
difference_geometries <- function(x, y) {
  validate_geometry(
    st_difference(
      st_geometry(validate_geometry(x)),
      union_geometries(y)
    )
  )
}


# ------------------------------------------------------------------------------
# as_sf_geometry()
# ------------------------------------------------------------------------------
# Purpose:
#   Wrap an sfc geometry into an sf object with only a geometry column.
# ------------------------------------------------------------------------------
as_sf_geometry <- function(x) {
  st_sf(geometry = validate_geometry(x))
}


# ==============================================================================
# 4. File input helper
# ==============================================================================

# ------------------------------------------------------------------------------
# read_polygon_file()
# ------------------------------------------------------------------------------
# Purpose:
#   Read an uploaded polygon file into an sf object.
#
# Supported formats:
#   - GeoJSON: .geojson
#   - JSON:    .json
#   - GeoPackage: .gpkg
#   - Zipped shapefile: .zip
#
# Important:
#   For zipped shapefiles, the ZIP must contain at least one .shp file together
#   with its supporting files, such as .dbf, .shx, and optionally .prj.
# ------------------------------------------------------------------------------
read_polygon_file <- function(file) {
  req(file)
  
  ext <- tolower(file_ext(file$name))
  
  if (ext %in% c("geojson", "json", "gpkg")) {
    
    poly <- st_read(file$datapath, quiet = TRUE)
    
  } else if (ext == "zip") {
    
    tmpdir <- tempfile("shp_")
    dir.create(tmpdir)
    
    unzip(file$datapath, exdir = tmpdir)
    
    shp_file <- list.files(
      tmpdir,
      pattern = "\\.shp$",
      full.names = TRUE,
      recursive = TRUE
    )
    
    if (length(shp_file) == 0) {
      stop("No .shp file was found inside the uploaded ZIP file.")
    }
    
    poly <- st_read(shp_file[1], quiet = TRUE)
    
  } else {
    
    stop(
      "Unsupported file type. Please upload GeoJSON, JSON, GPKG, ",
      "or a zipped shapefile."
    )
  }
  
  validate_geometry(poly)
}


# ==============================================================================
# 5. Leaflet drawing helper
# ==============================================================================

# ------------------------------------------------------------------------------
# leaflet_feature_to_sf()
# ------------------------------------------------------------------------------
# Purpose:
#   Convert a polygon drawn with leaflet.extras into an sf object.
#
# Supported geometries:
#   - Polygon
#   - MultiPolygon
#
# Limitations:
#   - Holes/interior rings are not fully preserved in this simplified converter.
#     If holes are needed in future versions, this function should be extended to
#     parse all coordinate rings, not only the first ring.
# ------------------------------------------------------------------------------
leaflet_feature_to_sf <- function(feature) {
  req(feature)
  
  geom_type <- feature$geometry$type
  
  if (geom_type == "Polygon") {
    
    coords <- feature$geometry$coordinates[[1]]
    
    lng <- vapply(coords, function(x) x[[1]], numeric(1))
    lat <- vapply(coords, function(x) x[[2]], numeric(1))
    
    geom <- st_polygon(list(cbind(lng, lat)))
    
  } else if (geom_type == "MultiPolygon") {
    
    polygons <- lapply(feature$geometry$coordinates, function(poly_coords) {
      coords <- poly_coords[[1]]
      
      lng <- vapply(coords, function(x) x[[1]], numeric(1))
      lat <- vapply(coords, function(x) x[[2]], numeric(1))
      
      st_polygon(list(cbind(lng, lat)))
    })
    
    geom <- st_multipolygon(polygons)
    
  } else {
    
    stop("Only Polygon and MultiPolygon geometries are supported.")
  }
  
  st_sf(geometry = st_sfc(geom, crs = CRS_LONG_LAT))
}


# ==============================================================================
# 6. Barrier construction
# ==============================================================================

# ------------------------------------------------------------------------------
# build_default_world_barrier()
# ------------------------------------------------------------------------------
# Purpose:
#   Build a default global land barrier from Natural Earth country polygons.
#
# Arguments:
#   barrier_sresol:
#     Resolution parameter used to define an approximate scale for filtering
#     small isolated polygons.
#
#   dist_simplify_km:
#     Simplification tolerance in kilometres. Higher values produce simpler
#     polygons and faster operations, but may remove narrow coastal details.
#
#   remove_small_isolated:
#     If TRUE, removes polygons that are both small and spatially isolated.
#
# Returns:
#   A valid sf object in EPSG:4326.
# ------------------------------------------------------------------------------
build_default_world_barrier <- function(
    barrier_sresol = DEFAULT_BARRIER_SRESOL,
    dist_simplify_km = 20,
    remove_small_isolated = TRUE
) {
  validate(
    need(barrier_sresol > 0, "Barrier resolution must be positive."),
    need(dist_simplify_km >= 0, "Simplification tolerance cannot be negative.")
  )
  
  max_edge_km <- EARTH_CIRCUMFERENCE_KM * 0.2 / barrier_sresol
  buffer_km <- max_edge_km / 2
  
  # Load global country polygons.
  world_ll <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf"
  )
  
  # Work in Mollweide kilometres for distance-like operations.
  world_moll <- transform_to_mollweide(st_geometry(world_ll))
  
  if (remove_small_isolated) {
    
    # Polygon area is used as an approximate size filter.
    world_area <- st_area(world_moll)
    
    # Neighbourhood graph: polygons with no neighbours are isolated.
    nb <- spdep::poly2nb(pl = world_moll, snap = buffer_km)
    number_of_neighbours <- spdep::card(nb)
    
    # Classify polygons by size and connectivity.
    polygon_class <- paste(
      ifelse(as.numeric(sqrt(world_area)) > max_edge_km, "big", "small"),
      ifelse(number_of_neighbours > 0, "connected", "isolated")
    )
    
    barrier_moll <- world_moll[polygon_class != "small isolated"]
    
  } else {
    
    barrier_moll <- world_moll
  }
  
  # Dissolve all polygons into one barrier object.
  barrier_moll <- union_geometries(barrier_moll)
  
  # Simplify geometry to improve performance.
  # preserveTopology = FALSE is faster, but may occasionally create geometry
  # artefacts; st_make_valid() is applied afterwards.
  barrier_moll <- st_simplify(
    barrier_moll,
    dTolerance = dist_simplify_km,
    preserveTopology = FALSE
  )
  
  transform_to_longlat(as_sf_geometry(barrier_moll))
}


# ==============================================================================
# 7. Spherical mesh construction
# ==============================================================================

# ------------------------------------------------------------------------------
# build_spherical_mesh()
# ------------------------------------------------------------------------------
# Purpose:
#   Build a spherical mesh over non-barrier locations and classify mesh triangles
#   whose centroids fall inside the barrier polygon.
#
# Arguments:
#   barrier_ll:
#     Barrier polygon in longitude/latitude coordinates.
#
#   mesh_sresol:
#     Spherical resolution used to generate candidate globe points.
#
#   globe_resol:
#     Resolution passed to fm_rcdt_2d_inla().
#
#   cutoff_multiplier:
#     Controls the minimum separation between mesh nodes.
#
# Returns:
#   A list containing:
#     - mesh: INLA/fmesher mesh object.
#     - tri_centers_ll: triangle centroids in EPSG:4326.
#     - tri_barrier: indices of barrier triangles.
#     - n_nodes: number of mesh nodes.
#     - n_triangles: number of mesh triangles.
#     - n_barrier_triangles: number of barrier triangles.
# ------------------------------------------------------------------------------
build_spherical_mesh <- function(
    barrier_ll,
    mesh_sresol = DEFAULT_MESH_SRESOL,
    globe_resol = 10,
    cutoff_multiplier = 1
) {
  validate(
    need(!is.null(barrier_ll), "A barrier polygon is required."),
    need(mesh_sresol > 0, "Mesh spherical resolution must be positive."),
    need(globe_resol > 0, "Globe resolution must be positive."),
    need(cutoff_multiplier > 0, "Cutoff multiplier must be positive.")
  )
  
  barrier_moll <- transform_to_mollweide(barrier_ll)
  
  buffer_km <- (EARTH_CIRCUMFERENCE_KM * 0.2 / mesh_sresol) / 2
  
  # Generate approximately regular points on the sphere.
  regular_globe_points <- fmesher_globe_points(globe = mesh_sresol)
  
  # Convert 3D spherical coordinates to sf points and transform them to the
  # Mollweide CRS so that they can be tested against the barrier polygon.
  initial_points_moll <- fm_transform(
    st_as_sf(
      as.data.frame(regular_globe_points),
      coords = 1:3,
      crs = fm_crs("sphere")
    ),
    fm_crs(barrier_moll)
  )
  
  # Remove candidate points that fall inside a slightly eroded land barrier.
  # The negative buffer helps avoid placing mesh nodes very close to land.
  barrier_for_filtering <- st_buffer(
    st_union(barrier_moll),
    dist = -buffer_km / 2
  )
  
  ocean_idx <- which(
    lengths(st_intersects(initial_points_moll, barrier_for_filtering)) == 0
  )
  
  if (length(ocean_idx) == 0) {
    stop("No ocean candidate points were found. Check the barrier geometry.")
  }
  
  # Build spherical triangulation over the retained ocean points.
  smesh <- fm_rcdt_2d_inla(
    loc = regular_globe_points[ocean_idx, ],
    globe = globe_resol,
    cutoff = cutoff_multiplier / mesh_sresol,
    crs = fm_crs("sphere")
  )
  
  # Identify triangles whose centroids are inside the barrier.
  tri_barrier <- unlist(
    fm_contains(
      x = barrier_moll,
      y = smesh,
      type = "centroid"
    )
  )
  
  # Compute triangle centroids in 3D Cartesian coordinates.
  tri_centers_xyz0 <- cbind(
    smesh$loc[smesh$graph$tv[, 1], 1:3] +
      smesh$loc[smesh$graph$tv[, 2], 1:3] +
      smesh$loc[smesh$graph$tv[, 3], 1:3]
  ) / 3
  
  # Normalize centroids back onto the unit sphere.
  tri_centers_xyz <- tri_centers_xyz0 / sqrt(rowSums(tri_centers_xyz0^2))
  
  # Convert triangle centroids to longitude/latitude for Leaflet display.
  tri_centers_ll <- fm_transform(
    st_as_sf(
      as.data.frame(tri_centers_xyz),
      coords = 1:3,
      crs = fm_crs(smesh)
    ),
    CRS_LONG_LAT
  )
  
  tri_centers_ll$inside_barrier <- FALSE
  tri_centers_ll$inside_barrier[tri_barrier] <- TRUE
  
  list(
    mesh = smesh,
    tri_centers_ll = tri_centers_ll,
    tri_barrier = tri_barrier,
    n_nodes = smesh$n,
    n_triangles = nrow(smesh$graph$tv),
    n_barrier_triangles = length(tri_barrier)
  )
}


# ==============================================================================
# 8. Local correlation diagnostic
# ==============================================================================

# ------------------------------------------------------------------------------
# localCorrel2D()
# ------------------------------------------------------------------------------
# Purpose:
#   Compute the local correlation between a selected location and all mesh nodes.
#
# Arguments:
#   locs_ll:
#     Matrix of selected locations in longitude/latitude coordinates.
#
#   mesh_locs_ll:
#     Mesh node coordinates in longitude/latitude.
#
#   Q:
#     Sparse precision matrix from the barrier SPDE model.
#
# Method:
#   1. Find the closest mesh node to each selected location using s2 distances.
#   2. Solve Q^{-1} b, where b selects the nearest mesh node.
#   3. Standardize covariance values into correlations.
# ------------------------------------------------------------------------------
localCorrel2D <- function(locs_ll, mesh_locs_ll, Q) {
  nl <- nrow(locs_ll)
  
  mesh_s2 <- s2::s2_lnglat(mesh_locs_ll[, 1], mesh_locs_ll[, 2])
  
  nearest_node_idx <- sapply(seq_len(nl), function(i) {
    loc_s2 <- s2::s2_lnglat(locs_ll[i, 1], locs_ll[i, 2])
    
    which.min(
      s2::s2_distance(
        mesh_s2,
        loc_s2,
        radius = s2::s2_earth_radius_meters()
      )
    )
  })
  
  # Indicator matrix selecting the nearest node for each clicked location.
  b <- matrix(0, nrow(Q), nl)
  
  for (i in seq_len(nl)) {
    b[nearest_node_idx[i], i] <- 1
  }
  
  # Solve Q x = b.
  covariance_to_selected_node <- inla.qsolve(Q, b)
  
  # Marginal standard deviations from Q^{-1}.
  marginal_sd <- sqrt(diag(inla.qinv(Q)))
  
  # Convert covariance to correlation.
  for (i in seq_len(nl)) {
    covariance_to_selected_node[, i] <-
      covariance_to_selected_node[, i] /
      (marginal_sd * marginal_sd[nearest_node_idx[i]])
  }
  
  drop(covariance_to_selected_node)
}


# ------------------------------------------------------------------------------
# build_correlation_grid()
# ------------------------------------------------------------------------------
# Purpose:
#   Project mesh-node correlations onto a regular longitude/latitude grid for
#   display on Leaflet as coloured rectangles.
#
# Arguments:
#   smesh:
#     Spherical mesh.
#
#   mcorrels:
#     Correlation values defined on mesh nodes.
#
#   grid_resol_deg:
#     Resolution of the display grid in degrees.
#
#   threshold:
#     Minimum correlation value to plot.
#
#   corr_id:
#     Identifier for the selected point. Allows several correlation surfaces to
#     remain visible simultaneously.
# ------------------------------------------------------------------------------
build_correlation_grid <- function(
    smesh,
    mcorrels,
    grid_resol_deg = 2,
    threshold = 0.1,
    corr_id = 1
) {
  validate(
    need(grid_resol_deg > 0, "Grid resolution must be positive."),
    need(threshold >= 0 && threshold <= 1, "Threshold must be between 0 and 1.")
  )
  
  lon_seq <- seq(-180, 180, by = grid_resol_deg)
  lat_seq <- seq(-85, 85, by = grid_resol_deg)
  
  grid_df <- expand.grid(lon = lon_seq, lat = lat_seq)
  
  # Convert longitude/latitude grid to spherical coordinates.
  grid_spherical <- inla.mesh.map(
    as.matrix(grid_df[, c("lon", "lat")]),
    projection = "longlat",
    inverse = TRUE
  )
  
  # Build projector from mesh nodes to grid points.
  grid_projector <- inla.mesh.projector(
    mesh = smesh,
    loc = grid_spherical
  )
  
  # Transform correlations before interpolation to keep values inside (-1, 1).
  field_transformed <- qlogis(
    0.5 + (0.5 - 1e-9) * as.numeric(mcorrels)
  )
  
  projected_correlations <- -1 + 2 * plogis(
    as.numeric(
      inla.mesh.project(
        projector = grid_projector,
        field = field_transformed
      )
    )
  )
  
  grid_df$correlation <- projected_correlations
  
  grid_df <- grid_df[is.finite(grid_df$correlation), ]
  grid_df <- grid_df[grid_df$correlation > threshold, ]
  
  if (nrow(grid_df) == 0) {
    stop("No correlations above the plotting threshold.")
  }
  
  # Rectangle bounds for Leaflet.
  grid_df$xmin <- grid_df$lon - grid_resol_deg / 2
  grid_df$xmax <- grid_df$lon + grid_resol_deg / 2
  grid_df$ymin <- grid_df$lat - grid_resol_deg / 2
  grid_df$ymax <- grid_df$lat + grid_resol_deg / 2
  
  grid_df$corr_id <- corr_id
  
  grid_df
}


# ==============================================================================
# 9. User interface
# ==============================================================================

ui <- navbarPage(
  
  title = HTML("<b>MOBIE</b>: Modifying Ocean Barriers for INLA on Earth"),
  
  tabPanel(
    "App",
    
    sidebarLayout(
      
      sidebarPanel(
        
        # ----------------------------------------------------------------------
        # 1. Barrier source
        # ----------------------------------------------------------------------
        h4("1. Barrier source"),
        
        fileInput(
          "file",
          "Upload polygon file",
          accept = c(".geojson", ".json", ".gpkg", ".zip")
        ),
        
        hr(),
        
        # ----------------------------------------------------------------------
        # 2. Barrier construction
        # ----------------------------------------------------------------------
        h4("2. Barrier construction"),
        
        numericInput(
          "barrier_sresol",
          "Barrier base spherical resolution",
          value = DEFAULT_BARRIER_SRESOL,
          min = 5,
          max = 100,
          step = 5
        ),
        
        numericInput(
          "simplify",
          "Simplification tolerance, km",
          value = 20,
          min = 0,
          step = 5
        ),
        
        checkboxInput(
          "remove_small_isolated",
          "Remove small isolated islands",
          TRUE
        ),
        
        actionButton(
          "build_default",
          "Build default world barrier",
          class = "btn-primary"
        ),
        
        hr(),
        
        # ----------------------------------------------------------------------
        # 3. Manual editing
        # ----------------------------------------------------------------------
        h4("3. Manual editing"),
        
        actionButton(
          "cut",
          "Cut drawn polygon",
          class = "btn-danger"
        ),
        
        actionButton(
          "add",
          "Add drawn polygon",
          class = "btn-success"
        ),
        
        hr(),
        
        # ----------------------------------------------------------------------
        # 4. Mesh construction
        # ----------------------------------------------------------------------
        h4("4. Mesh"),
        
        numericInput(
          "mesh_sresol",
          "Mesh base spherical resolution",
          value = DEFAULT_MESH_SRESOL,
          min = 5,
          max = 100,
          step = 5
        ),
        
        numericInput(
          "globe_resol",
          "Mesh globe resolution inside barriers",
          value = 10,
          min = 3,
          max = 40,
          step = 1
        ),
        
        numericInput(
          "cutoff_multiplier",
          "Mesh cutoff multiplier",
          value = 1,
          min = 0.1,
          max = 5,
          step = 0.1
        ),
        
        actionButton(
          "build_mesh",
          "Build spherical mesh",
          class = "btn-warning"
        ),
        
        hr(),
        
        # ----------------------------------------------------------------------
        # 5. Correlation diagnostic
        # ----------------------------------------------------------------------
        h4("5. Correlation diagnostic"),
        
        numericInput(
          "corr_range",
          "Correlation range",
          value = 0.3,
          min = 0.01,
          step = 0.05
        ),
        
        numericInput(
          "corr_sigma",
          "Correlation sigma",
          value = 1,
          min = 0.01,
          step = 0.1
        ),
        
        numericInput(
          "corr_grid_resol",
          "Correlation grid resolution, degrees",
          value = 2,
          min = 0.5,
          max = 10,
          step = 0.5
        ),
        
        numericInput(
          "corr_threshold",
          "Minimum correlation to plot",
          value = 0.1,
          min = 0,
          max = 1,
          step = 0.05
        ),
        
        actionButton(
          "clear_corr",
          "Clear correlations",
          class = "btn-secondary"
        ),
        
        hr(),
        
        # ----------------------------------------------------------------------
        # Downloads
        # ----------------------------------------------------------------------
        downloadButton(
          "download_barrier",
          "Download barrier GeoJSON"
        ),
        
        downloadButton(
          "download_mesh_rdata",
          "Download mesh .RData"
        ),
        
        downloadButton(
          "download_parameters_txt",
          "Download parameters .txt"
        ),
        
        hr(),
        
        # ----------------------------------------------------------------------
        # Status panel
        # ----------------------------------------------------------------------
        verbatimTextOutput("status")
      ),
      
      mainPanel(
        leafletOutput("map", height = "760px")
      )
    )
  ),
  
  tabPanel(
    "Documentation",
    
    fluidPage(
      
      h2("MOBIE documentation"),
      
      p(
        "MOBIE is an interactive tool for creating, editing, visualising, ",
        "and exporting barrier polygons and spherical meshes for INLA-based ",
        "global spatial models."
      ),
      
      h3("1. Barrier source"),
      h4("Upload polygon file"),
      p(
        "Upload an existing barrier polygon in GeoJSON, JSON, GeoPackage, ",
        "or zipped shapefile format."
      ),
      
      h3("2. Barrier construction"),
      h4("Barrier base spherical resolution"),
      p(
        "Controls the approximate resolution used when identifying small ",
        "isolated land polygons during construction of the default barrier."
      ),
      
      h4("Simplification tolerance, km"),
      p(
        "Controls how strongly the barrier geometry is simplified. Larger ",
        "values reduce geometry complexity and improve performance, but may ",
        "remove important narrow features."
      ),
      
      h4("Remove small isolated islands"),
      p(
        "If enabled, small isolated land polygons are removed before the ",
        "barrier is dissolved and simplified."
      ),
      
      h4("Build default world barrier"),
      p(
        "Creates a default global land barrier from Natural Earth country ",
        "polygons."
      ),
      
      h3("3. Manual editing"),
      h4("Draw polygon"),
      p(
        "Use the drawing tool on the map to create a polygon. The drawn ",
        "polygon can then be added to or subtracted from the current barrier."
      ),
      
      h4("Cut drawn polygon"),
      p(
        "Subtracts the drawn polygon from the current barrier."
      ),
      
      h4("Add drawn polygon"),
      p(
        "Merges the drawn polygon with the current barrier."
      ),
      
      h3("4. Mesh"),
      h4("Mesh base spherical resolution"),
      p(
        "Controls the density of candidate locations used to build the ",
        "spherical mesh."
      ),
      
      h4("Mesh globe resolution inside barriers"),
      p(
        "Controls the resolution parameter passed to the spherical mesh ",
        "constructor."
      ),
      
      h4("Mesh cutoff multiplier"),
      p(
        "Controls the minimum separation between retained mesh nodes."
      ),
      
      h4("Build spherical mesh"),
      p(
        "Builds the spherical mesh and identifies mesh triangles whose ",
        "centroids fall inside the barrier polygon."
      ),
      
      h3("5. Correlation diagnostic"),
      
      p(
        "This diagnostic visualizes the spatial correlation structure implied by ",
        "the barrier-SPDE model. After building the mesh, clicking on a location ",
        "computes the correlation between that location and all other locations ",
        "on the globe, accounting for the barrier geometry."
      ),
      
      h4("Correlation range"),
      
      p(
        "Controls the spatial scale of dependence in the barrier-SPDE model. ",
        "Larger values produce broader and smoother correlation patterns, whereas ",
        "smaller values generate more localized correlation structures. ",
        "The range parameter determines how rapidly correlation decays with ",
        "distance."
      ),
      
      h4("Correlation sigma"),
      
      p(
        "Controls the marginal standard deviation of the underlying spatial field. ",
        "Although correlation values are standardized, sigma influences the ",
        "construction of the precision matrix used to compute the diagnostic."
      ),
      
      h4("Correlation grid resolution"),
      
      p(
        "Defines the spatial resolution of the projected correlation map. ",
        "Smaller values produce finer maps but require more computation and memory. ",
        "Larger values generate coarser but faster visualizations."
      ),
      
      h4("Minimum correlation to plot"),
      
      p(
        "Only correlations above this threshold are displayed. Increasing the ",
        "threshold highlights the strongest dependencies and reduces visual clutter."
      ),
      
      h4("Click on the map"),
      
      p(
        "Each click selects a reference location and computes the correlation ",
        "between that location and all mesh nodes. The resulting correlation ",
        "surface is projected onto a regular longitude-latitude grid and added ",
        "to the map. Multiple correlation surfaces can be displayed simultaneously."
      ),
      
      h4("Clear correlations"),
      
      p(
        "Removes all correlation surfaces and all selected reference locations."
      ),
      
      h3("Downloads"),
      h4("Download barrier GeoJSON"),
      p(
        "Exports the current edited barrier polygon as a GeoJSON file."
      ),
      
      h4("Download mesh .RData"),
      p(
        "Exports the spherical mesh object, called smesh, and the vector of ",
        "barrier-triangle indices, called tri_barrier."
      )
    )
  )
)


# ==============================================================================
# 10. Server
# ==============================================================================

server <- function(input, output, session) {
  
  # ---------------------------------------------------------------------------
  # Reactive state
  # ---------------------------------------------------------------------------
  
  # Current barrier polygon as an sf object in EPSG:4326.
  barrier <- reactiveVal(NULL)
  
  # Most recent polygon drawn by the user on the Leaflet map.
  drawn_polygon <- reactiveVal(NULL)
  
  # Mesh construction result returned by build_spherical_mesh().
  mesh_data <- reactiveVal(NULL)
  
  # sf object containing clicked points used for correlation diagnostics.
  corr_points <- reactiveVal(NULL)
  
  # data.frame containing projected correlation grids for all clicked points.
  corr_grids <- reactiveVal(NULL)
  
  # Integer counter used to assign a unique ID to each clicked correlation point.
  corr_counter <- reactiveVal(0)
  
  
  # ---------------------------------------------------------------------------
  # Initial map
  # ---------------------------------------------------------------------------
  
  output$map <- renderLeaflet({
    
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) %>%
      
      addProviderTiles("CartoDB.Positron") %>%
      
      setView(lng = 0, lat = 20, zoom = 2) %>%
      
      addDrawToolbar(
        targetGroup = "drawn",
        
        polygonOptions = drawPolygonOptions(
          shapeOptions = drawShapeOptions(color = "#E67E22")
        ),
        
        editOptions = FALSE,
        polylineOptions = FALSE,
        circleOptions = FALSE,
        rectangleOptions = FALSE,
        markerOptions = FALSE,
        circleMarkerOptions = FALSE
      ) %>%
      
      addLayersControl(
        overlayGroups = c(
          "barrier",
          "drawn",
          "triangle centroids",
          "barrier triangle centroids",
          "correlation grid",
          "correlation point"
        ),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  
  # ---------------------------------------------------------------------------
  # refresh_map()
  # ---------------------------------------------------------------------------
  # Purpose:
  #   Redraw all map layers from the current reactive state.
  #
  # Why:
  #   Several actions modify the barrier, mesh, drawn polygon, or correlations.
  #   Centralising the map redraw avoids duplicated map-rendering logic.
  # ---------------------------------------------------------------------------
  
  refresh_map <- function() {
    
    proxy <- leafletProxy("map") %>%
      clearGroup("barrier") %>%
      clearGroup("drawn") %>%
      clearGroup("triangle centroids") %>%
      clearGroup("barrier triangle centroids") %>%
      clearGroup("correlation grid") %>%
      clearGroup("correlation point")
    
    # Draw current barrier.
    if (!is.null(barrier())) {
      proxy <- proxy %>%
        addPolygons(
          data = barrier(),
          group = "barrier",
          fillOpacity = 0.35,
          color = "#1f78b4",
          weight = 0.5,
          fillColor = "#1f78b4",
          label = "Barrier polygon"
        )
    }
    
    # Draw current user-drawn polygon.
    if (!is.null(drawn_polygon())) {
      proxy <- proxy %>%
        addPolygons(
          data = drawn_polygon(),
          group = "drawn",
          fillOpacity = 0.35,
          color = "#e67e22",
          weight = 2,
          fillColor = "#e67e22",
          label = "Drawn polygon"
        )
    }
    
    # Draw mesh triangle centroids, separated by ocean/barrier classification.
    if (!is.null(mesh_data())) {
      
      md <- mesh_data()
      
      tri <- md$tri_centers_ll
      tri_ocean <- tri[!tri$inside_barrier, ]
      tri_barrier <- tri[tri$inside_barrier, ]
      
      if (nrow(tri_ocean) > 0) {
        proxy <- proxy %>%
          addCircleMarkers(
            data = tri_ocean,
            group = "triangle centroids",
            radius = 2,
            stroke = FALSE,
            fillOpacity = 0.25,
            fillColor = "#4daf4a",
            label = "Ocean triangle centroid"
          )
      }
      
      if (nrow(tri_barrier) > 0) {
        proxy <- proxy %>%
          addCircleMarkers(
            data = tri_barrier,
            group = "barrier triangle centroids",
            radius = 2,
            stroke = FALSE,
            fillOpacity = 0.65,
            fillColor = "#000000",
            label = "Barrier triangle centroid"
          )
      }
    }
    
    # Draw projected correlation grids.
    if (!is.null(corr_grids())) {
      
      all_corr <- corr_grids()
      
      pal <- colorNumeric(
        palette = "YlOrRd",
        domain = all_corr$correlation
      )
      
      proxy <- proxy %>%
        addRectangles(
          data = all_corr,
          lng1 = ~xmin,
          lat1 = ~ymin,
          lng2 = ~xmax,
          lat2 = ~ymax,
          group = "correlation grid",
          fillColor = ~pal(correlation),
          fillOpacity = 0.55,
          stroke = FALSE,
          label = ~paste0(
            "Point ",
            corr_id,
            " | Correlation: ",
            round(correlation, 3)
          )
        )
    }
    
    # Draw clicked correlation points.
    if (!is.null(corr_points())) {
      proxy <- proxy %>%
        addCircleMarkers(
          data = corr_points(),
          group = "correlation point",
          radius = 7,
          color = "white",
          fillColor = "black",
          fillOpacity = 1,
          weight = 2,
          label = ~paste0("Correlation point ", corr_id)
        )
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # reset_correlations()
  # ---------------------------------------------------------------------------
  # Purpose:
  #   Clear all correlation diagnostic state.
  #
  # When used:
  #   - When the barrier changes.
  #   - When the mesh changes.
  #   - When the user clicks "Clear correlations".
  # ---------------------------------------------------------------------------
  
  reset_correlations <- function() {
    corr_points(NULL)
    corr_grids(NULL)
    corr_counter(0)
  }
  
  
  # ---------------------------------------------------------------------------
  # Upload polygon file
  # ---------------------------------------------------------------------------
  
  observeEvent(input$file, {
    
    req(input$file)
    
    tryCatch({
      
      poly <- read_polygon_file(input$file)
      
      # Dissolve uploaded polygons into one barrier geometry in Mollweide,
      # then convert back to longitude/latitude for map display and storage.
      poly <- as_sf_geometry(
        union_geometries(transform_to_mollweide(poly))
      )
      
      barrier(transform_to_longlat(poly))
      
      # Existing mesh and diagnostics are no longer valid after barrier changes.
      mesh_data(NULL)
      drawn_polygon(NULL)
      reset_correlations()
      
      refresh_map()
      
      showNotification(
        "Polygon uploaded successfully.",
        type = "message"
      )
      
    }, error = function(e) {
      
      showNotification(e$message, type = "error")
    })
  })
  
  
  # ---------------------------------------------------------------------------
  # Build default world barrier
  # ---------------------------------------------------------------------------
  
  observeEvent(input$build_default, {
    
    withProgress(message = "Building default world barrier...", value = 0.3, {
      
      tryCatch({
        
        poly <- build_default_world_barrier(
          barrier_sresol = input$barrier_sresol,
          dist_simplify_km = input$simplify,
          remove_small_isolated = input$remove_small_isolated
        )
        
        incProgress(0.7)
        
        barrier(poly)
        
        # Existing mesh and diagnostics are no longer valid after barrier changes.
        mesh_data(NULL)
        drawn_polygon(NULL)
        reset_correlations()
        
        refresh_map()
        
        showNotification(
          "Default barrier created.",
          type = "message"
        )
        
      }, error = function(e) {
        
        showNotification(e$message, type = "error")
      })
    })
  })
  
  
  # ---------------------------------------------------------------------------
  # Store newly drawn polygon
  # ---------------------------------------------------------------------------
  
  observeEvent(input$map_draw_new_feature, {
    
    tryCatch({
      
      drawn_polygon(
        leaflet_feature_to_sf(input$map_draw_new_feature)
      )
      
      refresh_map()
      
    }, error = function(e) {
      
      showNotification(e$message, type = "error")
    })
  })
  
  
  # ---------------------------------------------------------------------------
  # Cut drawn polygon from barrier
  # ---------------------------------------------------------------------------
  
  observeEvent(input$cut, {
    
    req(barrier(), drawn_polygon())
    
    tryCatch({
      
      result <- difference_geometries(
        transform_to_mollweide(barrier()),
        transform_to_mollweide(drawn_polygon())
      )
      
      barrier(transform_to_longlat(as_sf_geometry(result)))
      
      # Existing mesh and diagnostics are no longer valid after barrier changes.
      drawn_polygon(NULL)
      mesh_data(NULL)
      reset_correlations()
      
      refresh_map()
      
      showNotification(
        "Polygon cut from barrier.",
        type = "message"
      )
      
    }, error = function(e) {
      
      showNotification(e$message, type = "error")
    })
  })
  
  
  # ---------------------------------------------------------------------------
  # Add drawn polygon to barrier
  # ---------------------------------------------------------------------------
  
  observeEvent(input$add, {
    
    req(barrier(), drawn_polygon())
    
    tryCatch({
      
      result <- union_geometries(
        list(
          transform_to_mollweide(barrier()),
          transform_to_mollweide(drawn_polygon())
        )
      )
      
      barrier(transform_to_longlat(as_sf_geometry(result)))
      
      # Existing mesh and diagnostics are no longer valid after barrier changes.
      drawn_polygon(NULL)
      mesh_data(NULL)
      reset_correlations()
      
      refresh_map()
      
      showNotification(
        "Polygon added to barrier.",
        type = "message"
      )
      
    }, error = function(e) {
      
      showNotification(e$message, type = "error")
    })
  })
  
  
  # ---------------------------------------------------------------------------
  # Build spherical mesh
  # ---------------------------------------------------------------------------
  
  observeEvent(input$build_mesh, {
    
    req(barrier())
    
    withProgress(message = "Building spherical mesh...", value = 0.2, {
      
      tryCatch({
        
        md <- build_spherical_mesh(
          barrier_ll = barrier(),
          mesh_sresol = input$mesh_sresol,
          globe_resol = input$globe_resol,
          cutoff_multiplier = input$cutoff_multiplier
        )
        
        incProgress(0.8)
        
        mesh_data(md)
        
        # Correlations depend on the mesh and must be recomputed.
        reset_correlations()
        
        refresh_map()
        
        showNotification(
          "Mesh created and classified.",
          type = "message"
        )
        
      }, error = function(e) {
        
        showNotification(e$message, type = "error")
      })
    })
  })
  
  
  # ---------------------------------------------------------------------------
  # Compute local correlation from map click
  # ---------------------------------------------------------------------------
  
  observeEvent(input$map_click, {
    
    req(barrier(), mesh_data())
    
    id <- corr_counter() + 1
    corr_counter(id)
    
    click <- input$map_click
    
    locs_ll <- matrix(
      c(click$lng, click$lat),
      ncol = 2,
      byrow = TRUE
    )
    
    selected_point <- st_as_sf(
      data.frame(
        lon = click$lng,
        lat = click$lat,
        corr_id = id
      ),
      coords = c("lon", "lat"),
      crs = CRS_LONG_LAT
    )
    
    withProgress(
      message = "Computing and projecting local correlation...",
      value = 0.2,
      {
        
        tryCatch({
          
          md <- mesh_data()
          
          smesh <- md$mesh
          tri_barrier <- md$tri_barrier
          
          # Build barrier finite-element matrices.
          fem <- inla.barrier.fem(
            mesh = smesh,
            barrier.triangles = tri_barrier
          )
          
          # Build barrier SPDE precision matrix.
          Q <- inla.barrier.q(
            fem = fem,
            ranges = c(input$corr_range, input$corr_range * 0.2),
            sigma = input$corr_sigma
          )
          
          # Convert mesh nodes to longitude/latitude for nearest-node search.
          mesh_coords <- inla.mesh.map(
            smesh$loc,
            projection = "longlat",
            inverse = FALSE
          )
          
          mcorrels <- localCorrel2D(
            locs_ll = locs_ll,
            mesh_locs_ll = mesh_coords,
            Q = Q
          )
          
          new_grid <- build_correlation_grid(
            smesh = smesh,
            mcorrels = mcorrels,
            grid_resol_deg = input$corr_grid_resol,
            threshold = input$corr_threshold,
            corr_id = id
          )
          
          if (is.null(corr_grids())) {
            corr_grids(new_grid)
          } else {
            corr_grids(rbind(corr_grids(), new_grid))
          }
          
          if (is.null(corr_points())) {
            corr_points(selected_point)
          } else {
            corr_points(rbind(corr_points(), selected_point))
          }
          
          refresh_map()
          
          showNotification(
            paste("Correlation point", id, "added."),
            type = "message"
          )
          
        }, error = function(e) {
          
          # Roll back the ID counter if the diagnostic fails.
          corr_counter(id - 1)
          
          showNotification(e$message, type = "error")
        })
      }
    )
  })
  
  
  # ---------------------------------------------------------------------------
  # Clear all correlation diagnostics
  # ---------------------------------------------------------------------------
  
  observeEvent(input$clear_corr, {
    
    reset_correlations()
    refresh_map()
    
    showNotification(
      "Correlation diagnostics cleared.",
      type = "message"
    )
  })
  
  
  # ---------------------------------------------------------------------------
  # Status output
  # ---------------------------------------------------------------------------
  
  output$status <- renderPrint({
    
    cat("Barrier loaded:", !is.null(barrier()), "\n")
    
    if (!is.null(barrier())) {
      cat("Barrier features:", nrow(barrier()), "\n")
      cat("Barrier CRS:", st_crs(barrier())$input, "\n")
    }
    
    cat("\nCurrent parameters\n")
    cat("Barrier base spherical resolution:", input$barrier_sresol, "\n")
    cat("Simplification tolerance, km:", input$simplify, "\n")
    cat("Remove small isolated islands:", input$remove_small_isolated, "\n")
    cat("Mesh base spherical resolution:", input$mesh_sresol, "\n")
    cat("Mesh globe resolution:", input$globe_resol, "\n")
    cat("Mesh cutoff multiplier:", input$cutoff_multiplier, "\n")
    
    if (!is.null(drawn_polygon())) {
      cat("\nDrawn polygon: yes\n")
    } else {
      cat("\nDrawn polygon: no\n")
    }
    
    if (!is.null(mesh_data())) {
      md <- mesh_data()
      
      cat("\nMesh summary\n")
      cat("Nodes:", md$n_nodes, "\n")
      cat("Triangles:", md$n_triangles, "\n")
      cat("Barrier triangles:", md$n_barrier_triangles, "\n")
    } else {
      cat("\nMesh: not built\n")
    }
    
    cat("\nCorrelation points:", corr_counter(), "\n")
  })
  
  
  # ---------------------------------------------------------------------------
  # Download edited barrier as GeoJSON
  # ---------------------------------------------------------------------------
  
  output$download_barrier <- downloadHandler(
    
    filename = function() {
      paste0("barrier_polygon_", Sys.Date(), ".geojson")
    },
    
    content = function(file) {
      
      req(barrier())
      
      st_write(
        barrier(),
        file,
        driver = "GeoJSON",
        delete_dsn = TRUE,
        quiet = TRUE
      )
    }
  )
  
  
  # ---------------------------------------------------------------------------
  # Download mesh as RData
  # ---------------------------------------------------------------------------
  
  output$download_mesh_rdata <- downloadHandler(
    
    filename = function() {
      paste0("MOBIE_mesh_", Sys.Date(), ".RData")
    },
    
    content = function(file) {
      
      req(mesh_data())
      
      md <- mesh_data()
      
      smesh <- md$mesh
      tri_barrier <- md$tri_barrier
      
      save(smesh, tri_barrier, file = file)
    }
  )

# ---------------------------------------------------------------------------
# Download parameters txt
# ---------------------------------------------------------------------------

output$download_parameters_txt <- downloadHandler(
  
  filename = function() {
    paste0("MOBIE_parameters_", Sys.Date(), ".txt")
  },
  
  content = function(file) {
    
    con <- file(file, open = "wt")
    on.exit(close(con), add = TRUE)
    
    writeLines("MOBIE: Modifying Ocean Barriers for INLA on Earth", con)
    writeLines("==================================================", con)
    writeLines("", con)
    
    writeLines("Export information", con)
    writeLines("------------------", con)
    writeLines(paste("Export date:", Sys.Date()), con)
    writeLines(paste("Export time:", format(Sys.time(), "%H:%M:%S")), con)
    writeLines("", con)
    
    writeLines("Barrier construction parameters", con)
    writeLines("-------------------------------", con)
    writeLines(
      paste(
        "Barrier base spherical resolution:",
        input$barrier_sresol
      ),
      con
    )
    
    writeLines(
      paste(
        "Simplification tolerance, km:",
        input$simplify
      ),
      con
    )
    
    writeLines(
      paste(
        "Remove small isolated islands:",
        input$remove_small_isolated
      ),
      con
    )
    
    writeLines("", con)
    
    writeLines("Mesh construction parameters", con)
    writeLines("----------------------------", con)
    
    writeLines(
      paste(
        "Mesh base spherical resolution:",
        input$mesh_sresol
      ),
      con
    )
    
    writeLines(
      paste(
        "Mesh globe resolution inside barriers:",
        input$globe_resol
      ),
      con
    )
    
    writeLines(
      paste(
        "Mesh cutoff multiplier:",
        input$cutoff_multiplier
      ),
      con
    )
    
    writeLines("", con)
    
    writeLines("Correlation diagnostic parameters", con)
    writeLines("---------------------------------", con)
    
    writeLines(
      paste(
        "Correlation range:",
        input$corr_range
      ),
      con
    )
    
    writeLines(
      paste(
        "Correlation sigma:",
        input$corr_sigma
      ),
      con
    )
    
    writeLines(
      paste(
        "Correlation grid resolution, degrees:",
        input$corr_grid_resol
      ),
      con
    )
    
    writeLines(
      paste(
        "Minimum correlation to plot:",
        input$corr_threshold
      ),
      con
    )
    
    writeLines("", con)
    
    writeLines("Barrier status", con)
    writeLines("--------------", con)
    
    writeLines(
      paste(
        "Barrier loaded:",
        !is.null(barrier())
      ),
      con
    )
    
    if (!is.null(barrier())) {
      
      writeLines(
        paste(
          "Barrier features:",
          nrow(barrier())
        ),
        con
      )
      
      writeLines(
        paste(
          "Barrier CRS:",
          st_crs(barrier())$input
        ),
        con
      )
    }
    
    writeLines("", con)
    
    writeLines("Mesh status", con)
    writeLines("-----------", con)
    
    if (!is.null(mesh_data())) {
      
      md <- mesh_data()
      
      writeLines("Mesh built: TRUE", con)
      
      writeLines(
        paste(
          "Number of mesh nodes:",
          md$n_nodes
        ),
        con
      )
      
      writeLines(
        paste(
          "Number of mesh triangles:",
          md$n_triangles
        ),
        con
      )
      
      writeLines(
        paste(
          "Number of barrier triangles:",
          md$n_barrier_triangles
        ),
        con
      )
      
    } else {
      
      writeLines("Mesh built: FALSE", con)
      
    }
    
    writeLines("", con)
    
    writeLines("Correlation status", con)
    writeLines("------------------", con)
    
    writeLines(
      paste(
        "Number of correlation points:",
        corr_counter()
      ),
      con
    )
    
    writeLines("", con)
    
    writeLines("Reproducibility information", con)
    writeLines("---------------------------", con)
    
    writeLines(
      paste(
        "R version:",
        R.version.string
      ),
      con
    )
    
    writeLines("", con)
    
    writeLines("Package versions", con)
    writeLines("----------------", con)
    
    packages <- c(
      "shiny",
      "leaflet",
      "leaflet.extras",
      "sf",
      "rnaturalearth",
      "spdep",
      "units",
      "s2",
      "INLA",
      "INLAspacetime",
      "inlabru"
    )
    
    for (pkg in packages) {
      
      version <- tryCatch(
        as.character(utils::packageVersion(pkg)),
        error = function(e) "not installed"
      )
      
      writeLines(
        paste(pkg, ":", version),
        con
      )
    }
  }
)

}

# ==============================================================================
# 11. Run app
# ==============================================================================

shinyApp(ui, server)