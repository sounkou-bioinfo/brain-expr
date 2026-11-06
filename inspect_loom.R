#!/usr/bin/env Rscript
# Usage : Rscript inspect_loom.R data/raw/dev_all.loom
# Affiche les attributs de cellules et leurs valeurs, pour ajuster config.R
# (noms d'attributs d'âge, de type cellulaire, de tissu).
suppressPackageStartupMessages(library(hdf5r))
f <- H5File$new(commandArgs(TRUE)[1], "r")
cat("matrix :", paste(f[["matrix"]]$dims, collapse = " x "), "\n")
for (a in names(f[["col_attrs"]])) {
  d <- f[[paste0("col_attrs/", a)]]
  if (length(d$dims) != 1) { cat(sprintf("%-22s [%s]\n", a, paste(d$dims, collapse = "x"))); next }
  v <- d[]
  u <- unique(v)
  cat(sprintf("%-22s %7d valeurs distinctes : %s\n", a, length(u),
              paste(head(if (is.numeric(u)) sort(u) else names(sort(table(v), decreasing = TRUE)), 12),
                    collapse = " | ")))
}
f$close_all()
