#' Run MOBIE
#'
#' Starts the MOBIE Shiny app.
#'
#' @param ... Extra arguments passed to [shiny::runApp()].
#'
#' @export
run_mobie <- function(...) {
  app_dir <- system.file("app", package = "mobie")

  if (identical(app_dir, "")) {
    app_dir <- file.path(getwd(), "inst", "app")
  }

  if (!dir.exists(app_dir)) {
    stop("Could not find the app directory. Run this from the MOBIE project root or install the package first.", call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
