#!/usr/bin/env Rscript
# Usage :  Rscript run_pipeline.R [config.R] [dossier_sortie]
# Produit dans results/ :
#   expression_summary.csv.gz  gène x type cellulaire x stade (moyenne, % cellules)
#   coexpression.csv.gz        paires de gènes x type cellulaire x stade
#   cells_per_group.csv        effectifs
#   dotplot_<dataset>_<niveau>.pdf
# puis copie les deux tables dans app/data/ pour l'app Shiny.

args <- commandArgs(trailingOnly = TRUE)
cfg_file <- if (length(args) >= 1) args[1] else "config.R"
out_dir  <- if (length(args) >= 2) args[2] else "results"

here <- function(...) file.path(dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), ...)
source(here("R", "loom_expr.R"))
source(cfg_file)
suppressPackageStartupMessages(library(ggplot2))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summ <- list(); coex <- list(); groups <- list()
for (nm in names(datasets)) {
  cfg <- datasets[[nm]]
  if (!file.exists(cfg$file)) { message("!! ", nm, " ignoré : ", cfg$file, " absent"); next }
  ds <- extract_dataset(nm, cfg, genes)
  if (is.null(ds)) next
  levels_avail <- sub("^ct_", "", grep("^ct_", names(ds$cells), value = TRUE))
  for (lv in levels_avail) {
    summ[[paste(nm, lv)]] <- summarise_expression(ds, lv, min_cells)
    coex[[paste(nm, lv)]] <- summarise_coexpression(ds, lv, min_cells)
    col <- paste0("ct_", lv)
    groups[[paste(nm, lv)]] <- ds$cells[, .(n_cells = .N, ages = paste(sort(unique(age)), collapse = ";")),
                                        by = .(dataset, stage, celltype = get(col))][, level := lv]
  }
}
if (!length(summ)) stop("Aucun jeu de données traité : vérifiez les chemins dans ", cfg_file)

summ <- rbindlist(summ); coex <- rbindlist(coex); groups <- rbindlist(groups)
fwrite(summ, file.path(out_dir, "expression_summary.csv.gz"))
fwrite(coex, file.path(out_dir, "coexpression.csv.gz"))
fwrite(groups, file.path(out_dir, "cells_per_group.csv"))

# --- dotplots statiques (un par jeu de données x niveau) ---
for (key in unique(summ[, paste(dataset, level)])) {
  d <- summ[paste(dataset, level) == key]
  d[, gene := factor(gene, levels = intersect(genes, unique(gene)))]
  ncat <- uniqueN(d$celltype)
  p <- ggplot(d, aes(gene, celltype)) +
    geom_point(aes(size = pct_expressing, colour = mean_lognorm)) +
    facet_wrap(~stage, nrow = 1) +
    scale_size_area(max_size = 6, name = "% cellules") +
    scale_colour_viridis_c(name = "moy. log(CP10k+1)") +
    labs(title = key, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
  ggsave(file.path(out_dir, sprintf("dotplot_%s.pdf", gsub(" ", "_", key))), p,
         width = 3 + 2.2 * uniqueN(d$stage) * length(levels(d$gene)) / 5,
         height = 2 + 0.16 * ncat, limitsize = FALSE)
}

app_data <- here("app", "data")
dir.create(app_data, showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(out_dir, c("expression_summary.csv.gz", "coexpression.csv.gz", "cells_per_group.csv")),
          app_data, overwrite = TRUE)
message("OK -> ", normalizePath(out_dir), " (copié aussi dans app/data/)")

# provenance affichée dans l'app
writeLines(c(paste("Généré le", format(Sys.time(), "%Y-%m-%d %H:%M")),
             paste("Config :", normalizePath(cfg_file)),
             paste("Gènes :", paste(genes, collapse = ", ")),
             vapply(names(datasets), function(n) sprintf("%s : %s (stades %s, fenêtre ±%s j)", n,
               datasets[[n]]$file, paste(datasets[[n]]$stages, collapse = ", "),
               datasets[[n]]$stage_window %||% 0), "")),
           file.path(app_data, "provenance.txt"))
