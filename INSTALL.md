# Install and run

From the project folder:

```r
devtools::load_all()
run_mobie()
```

Direct run:

```r
shiny::runApp("inst/app")
```

If something fails, first check that the same package versions used by the original script are installed.
