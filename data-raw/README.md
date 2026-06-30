# Data files

Large datasets are not stored in this repository.

The case studies require external input files such as cleaned GBIF data, FishGlob data, and MOBIE exports. These files should be stored locally or shared through a data repository such as Zenodo, OSF or Figshare.

Suggested local use:

```text
data-raw/
├── FishGlob_public_clean.csv
├── dat_clean_marine.RData
├── barrier_polygon_YYYY-MM-DD.geojson
└── MOBIE_mesh_YYYY-MM-DD.RData
```

If files are stored elsewhere, update the paths in the case-study scripts.
