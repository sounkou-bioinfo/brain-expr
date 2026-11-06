#!/usr/bin/env Rscript
# Builds one gene-major store per atlas (format: app/store.R) from the
# downloaded .loom/.h5ad files. Every gene is kept; the human atlases keep a
# fixed-seed random sample of cells per stage x cell class.
#
# Usage: Rscript extract/extract.R <raw_dir> <out_dir> [atlas ...]
#   atlases: lamanno2021 zeisel2018 braun2023 velmeshev2023 (default: all)

suppressPackageStartupMessages(library(rhdf5))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "..", "app", "store.R"))

CELLS_PER_GROUP <- 200L
SEED <- 1L

h5_text <- function(fid, path) as.character(h5read(fid, path))

# AnnData categorical column: integer codes into a categories array.
h5ad_factor <- function(fid, path) {
  x <- h5read(fid, path)
  if (!is.list(x)) return(as.character(x))
  as.character(x$categories)[x$codes + 1L]
}

sample_groups <- function(group, cap, seed = SEED) {
  set.seed(seed)
  idx <- split(seq_along(group), group, drop = TRUE)
  sort(unlist(lapply(idx, function(i) if (length(i) > cap) i[sample.int(length(i), cap)] else i),
              use.names = FALSE))
}

# Writes the store from per-gene sparse columns produced by `gene_columns`,
# a function(callback) that calls callback(gene_index, cell0, count) for every
# gene in gene order.
write_store <- function(out, meta, cells, genes, gene_columns) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  con <- file(file.path(out, "matrix.bin"), "wb")
  offset <- numeric(length(genes)); nbytes <- integer(length(genes)); nnz <- integer(length(genes))
  pos <- 0
  gene_columns(function(g, cell0, count) {
    blk <- encode_gene(cell0, count)
    writeBin(blk, con)
    offset[g] <<- pos; nbytes[g] <<- length(blk); nnz[g] <<- length(cell0)
    pos <<- pos + length(blk)
  })
  close(con)
  write_tsv_gz(cells, file.path(out, "cells.tsv.gz"))
  write_tsv_gz(data.frame(gene = genes, offset = sprintf("%.0f", offset), nbytes = nbytes, nnz = nnz),
               file.path(out, "genes.tsv.gz"))
  meta$levels <- paste(meta$levels, collapse = ",")
  meta$n_cells <- nrow(cells)
  meta$n_genes <- length(genes)
  meta$built <- format(Sys.time(), "%Y-%m-%d")
  write.dcf(as.data.frame(meta, check.names = FALSE), file.path(out, "meta.dcf"))
  message(sprintf("  %s: %d cells, %d genes, %.0f MB", basename(out), nrow(cells), length(genes), pos / 2^20))
}

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# --- loom: genes x cells, chunked; read bands of genes for the kept cells -----

extract_loom <- function(file, keep, cells, meta, out, band = 256L) {
  fid <- H5Fopen(file, "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))
  genes <- h5_text(fid, "row_attrs/Gene")
  total <- numeric(length(keep))
  write_store(out, meta, cells, genes, function(emit) {
    for (s in seq(1L, length(genes), by = band)) {
      g <- s:min(length(genes), s + band - 1L)
      m <- h5read(fid, "matrix", index = list(NULL, g))[keep, , drop = FALSE]
      total <<- total + rowSums(m)
      nz <- which(m > 0, arr.ind = TRUE)
      by_gene <- split(seq_len(nrow(nz)), factor(nz[, 2], seq_along(g)))
      for (j in seq_along(g)) {
        r <- nz[by_gene[[j]], 1]
        emit(g[j], r - 1L, m[cbind(r, j)])
      }
    }
  })
  cells$total <- total
  write_tsv_gz(cells, file.path(out, "cells.tsv.gz"))
}

# --- h5ad: cells x genes CSR; read the kept rows, then regroup by gene --------

read_csr_rows <- function(fid, xgrp, rows, batch = 2000L) {
  indptr <- h5read(fid, paste0(xgrp, "/indptr"), bit64conversion = "double")
  d <- H5Dopen(fid, paste0(xgrp, "/data")); di <- H5Dopen(fid, paste0(xgrp, "/indices"))
  on.exit({ H5Dclose(d); H5Dclose(di) })
  parts <- lapply(split(seq_along(rows), ceiling(seq_along(rows) / batch)), function(k) {
    r <- rows[k]
    n <- indptr[r + 1] - indptr[r]
    sp <- H5Dget_space(d); spi <- H5Dget_space(di)
    H5Sselect_none(sp); H5Sselect_none(spi)
    for (i in which(n > 0)) {
      H5Sselect_hyperslab(sp, start = indptr[r[i]] + 1, count = n[i], op = "H5S_SELECT_OR")
      H5Sselect_hyperslab(spi, start = indptr[r[i]] + 1, count = n[i], op = "H5S_SELECT_OR")
    }
    ms <- H5Screate_simple(sum(n))
    out <- list(cell = rep.int(k, n),
                gene = as.integer(H5Dread(di, h5spaceFile = spi, h5spaceMem = ms)) + 1L,
                count = as.integer(round(H5Dread(d, h5spaceFile = sp, h5spaceMem = ms))))
    H5Sclose(ms); H5Sclose(sp); H5Sclose(spi)
    out
  })
  list(cell = unlist(lapply(parts, `[[`, "cell")),
       gene = unlist(lapply(parts, `[[`, "gene")),
       count = unlist(lapply(parts, `[[`, "count")))
}

extract_h5ad <- function(file, xgrp, rows, genes, cells, meta, out) {
  fid <- H5Fopen(file, "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))
  t <- read_csr_rows(fid, xgrp, rows)
  total <- rowsum(as.numeric(t$count), t$cell)
  cells$total <- 0
  cells$total[as.integer(rownames(total))] <- total[, 1]
  o <- order(t$gene, t$cell)
  t <- lapply(t, `[`, o)
  ends <- cumsum(tabulate(t$gene, length(genes)))
  write_store(out, meta, cells, genes, function(emit) {
    start <- 1L
    for (g in seq_along(genes)) {
      k <- if (ends[g] >= start) start:ends[g] else integer()
      emit(g, t$cell[k] - 1L, t$count[k])
      start <- ends[g] + 1L
    }
  })
}

# --- atlases -----------------------------------------------------------------

atlas <- list(
  lamanno2021 = function(raw, out) {
    f <- file.path(raw, "dev_all.loom")
    fid <- H5Fopen(f, "H5F_ACC_RDONLY")
    age <- h5_text(fid, "col_attrs/Age"); cls <- h5_text(fid, "col_attrs/Class")
    sub <- h5_text(fid, "col_attrs/Subclass"); tis <- h5_text(fid, "col_attrs/Tissue")
    H5Fclose(fid)
    keep <- which(age %in% c("e12.5", "e13.5") & cls != "Bad cells")
    stage <- toupper(age[keep])
    cells <- data.frame(stage = stage, stage_rank = as.numeric(sub("E", "", stage)),
                        Class = cls[keep], Subclass = sub[keep], Tissue = tis[keep])
    extract_loom(f, keep, cells, out = out, meta = list(
      id = "lamanno2021", species = "mouse",
      label = "Souris embryonnaire E12.5-E13.5 (La Manno 2021)",
      citation = "La Manno et al. 2021, Nature 596:92-96. doi:10.1038/s41586-021-03775-x",
      levels = c("Class", "Subclass", "Tissue"),
      sampling = "toutes les cellules E12.5 et E13.5 (sauf Class = Bad cells)"))
  },

  zeisel2018 = function(raw, out) {
    f <- file.path(raw, "l5_all_rev1.loom")
    fid <- H5Fopen(f, "H5F_ACC_RDONLY")
    age <- h5_text(fid, "col_attrs/Age"); tis <- h5_text(fid, "col_attrs/Tissue")
    r3 <- h5_text(fid, "col_attrs/TaxonomyRank3"); r4 <- h5_text(fid, "col_attrs/TaxonomyRank4")
    H5Fclose(fid)
    in_window <- vapply(strsplit(age, ",\\s*"), function(a) all(a %in% sprintf("p%d", 19:23)), NA)
    keep <- which(in_window & !tis %in% c("DRG", "ENS", "Sympath", "SC"))
    cells <- data.frame(stage = "P19-P23", stage_rank = 21,
                        TaxonomyRank3 = r3[keep], TaxonomyRank4 = r4[keep], Tissue = tis[keep])
    extract_loom(f, keep, cells, out = out, meta = list(
      id = "zeisel2018", species = "mouse",
      label = "Souris P19-P23, cerveau (Zeisel 2018)",
      citation = "Zeisel et al. 2018, Cell 174:999-1014. doi:10.1016/j.cell.2018.06.021",
      levels = c("TaxonomyRank3", "TaxonomyRank4", "Tissue"),
      sampling = "toutes les cellules d'animaux P19-P23, sans moelle (SC), DRG, ENS ni ganglions sympathiques"))
  },

  braun2023 = function(raw, out) {
    f <- file.path(raw, "human_dev.h5ad")
    fid <- H5Fopen(f, "H5F_ACC_RDONLY")
    age <- round(as.numeric(h5ad_factor(fid, "obs/Age")), 1)
    cls <- h5ad_factor(fid, "obs/CellClass"); reg <- h5ad_factor(fid, "obs/Region")
    genes <- h5_text(fid, "var/Gene")
    H5Fclose(fid)
    rows <- sample_groups(paste(age, cls), CELLS_PER_GROUP)
    cells <- data.frame(stage = sprintf("%.1f pcw", age[rows]), stage_rank = age[rows],
                        CellClass = cls[rows], Region = reg[rows])
    extract_h5ad(f, "X", rows, genes, cells, out = out, meta = list(
      id = "braun2023", species = "human",
      label = "Humain 5-14 semaines post-conception, cerveau entier (Braun 2023)",
      citation = "Braun et al. 2023, Science 382:eadf1226. doi:10.1126/science.adf1226",
      levels = c("CellClass", "Region"),
      sampling = sprintf("au plus %d cellules par stade x CellClass, tirage aléatoire (graine %d)",
                         CELLS_PER_GROUP, SEED)))
  },

  velmeshev2023 = function(raw, out) {
    f <- file.path(raw, "velmeshev2023.h5ad")
    m <- read.delim(file.path(raw, "velmeshev_meta.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
    fid <- H5Fopen(f, "H5F_ACC_RDONLY")
    id <- h5_text(fid, "obs/_index"); genes <- h5ad_factor(fid, "raw/var/feature_name")
    H5Fclose(fid)
    m <- m[match(id, m$Cell_ID), ]
    stopifnot(!anyNA(m$Cell_ID))
    rank <- tapply(m[["Age_(days)"]], m$Age_Range, median)
    rows <- sample_groups(paste(m$Age_Range, m$Lineage), CELLS_PER_GROUP)
    cells <- data.frame(stage = m$Age_Range[rows], stage_rank = as.numeric(rank[m$Age_Range[rows]]),
                        Lineage = m$Lineage[rows], Region = m$Region_Broad[rows],
                        Study = m$Dataset[rows], age_days = m[["Age_(days)"]][rows])
    extract_h5ad(f, "raw/X", rows, genes, cells, out = out, meta = list(
      id = "velmeshev2023", species = "human",
      label = "Humain 2e trimestre - adulte, cortex (Velmeshev 2023)",
      citation = "Velmeshev et al. 2023, Science 382:eadf0834. doi:10.1126/science.adf0834",
      levels = c("Lineage", "Region", "Study"),
      sampling = sprintf("au plus %d cellules par tranche d'âge x Lineage, tirage aléatoire (graine %d)",
                         CELLS_PER_GROUP, SEED)))
  }
)

args <- commandArgs(TRUE)
raw <- args[1]; out <- args[2]
which_atlas <- if (length(args) > 2) args[-(1:2)] else names(atlas)
for (a in which_atlas) {
  message("== ", a)
  t0 <- Sys.time()
  atlas[[a]](raw, file.path(out, a))
  message("  ", format(Sys.time() - t0))
}
