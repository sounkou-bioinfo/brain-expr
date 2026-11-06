#!/usr/bin/env Rscript
# Builds the desktop app (Electron + webR) and the static web site from app/
# and the atlas stores.
#
# Usage: Rscript build_desktop.R <stores_dir> [win|mac|linux]
#        Rscript build_desktop.R <stores_dir> site     # static site only, in build/site
#
# The app runs R compiled to WebAssembly (shinylive), so the user installs
# nothing. The stores are served next to the app under data/ and read gene by
# gene with HTTP range requests.

args <- commandArgs(trailingOnly = TRUE)
stores <- normalizePath(args[1], mustWork = TRUE)
target <- if (length(args) > 1) args[2] else NULL
atlases <- c("lamanno2021", "zeisel2018", "braun2023", "velmeshev2023")

site <- file.path("build", "site")
unlink(site, recursive = TRUE)
shinylive::export("app", site, quiet = TRUE)
for (a in atlases) {
  dir.create(file.path(site, "data", a), recursive = TRUE)
  stopifnot(file.copy(list.files(file.path(stores, a), full.names = TRUE), file.path(site, "data", a)))
}
message("Site: ", site)
if (identical(target, "site")) quit(save = "no")

shinyelectron::build_electron_app(
  app_dir = site,
  output_dir = file.path("build", "electron-app"),
  app_name = "Brain Gene Expression",
  runtime_strategy = "shinylive",
  platform = target,
  overwrite = TRUE
)
message("Executables: build/electron-app/dist/")
