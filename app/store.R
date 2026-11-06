# Gene-major expression store, shared by the extraction scripts and the app.
#
# One directory per atlas:
#   meta.dcf      atlas label, species, citation, annotation levels, sampling
#   cells.tsv.gz  one row per kept cell: stage, stage_rank, annotation levels, total UMI
#   genes.tsv.gz  gene, offset, nbytes, nnz: where each gene's block sits in matrix.bin
#   matrix.bin    concatenated zlib blocks, one per gene
#
# A gene block holds the gene's non-zero cells: the 0-based cell indices
# delta-encoded as little-endian int32, then the counts as little-endian uint16.
# Base R reads it with readBin() and memDecompress(), so the same code runs in
# native R and in webR (shinylive, shinyelectron).

encode_gene <- function(cell0, count) {
  if (!length(cell0)) return(raw())
  memCompress(c(writeBin(as.integer(diff(c(0L, cell0))), raw(), size = 4, endian = "little"),
                writeBin(as.integer(pmin(count, 65535)), raw(), size = 2, endian = "little")),
              "gzip")
}

decode_gene <- function(block) {
  x <- memDecompress(block, "gzip")
  n <- length(x) %/% 6L
  list(cell = cumsum(readBin(x[seq_len(4L * n)], "integer", n, size = 4, endian = "little")) + 1L,
       count = readBin(x[4L * n + seq_len(2L * n)], "integer", n, size = 2,
                       signed = FALSE, endian = "little"))
}

# A store source is a local directory or a URL prefix ending in "/".
is_url <- function(src) grepl("^https?://", src)

store_file <- function(src, name) {
  if (!is_url(src)) return(file.path(src, name))
  tf <- tempfile(fileext = paste0("-", name))
  utils::download.file(paste0(src, name), tf, mode = "wb", quiet = TRUE)
  tf
}

store_range <- function(src, offset, nbytes) {
  if (!is_url(src)) {
    con <- file(file.path(src, "matrix.bin"), "rb")
    on.exit(close(con))
    seek(con, offset)
    return(readBin(con, "raw", nbytes))
  }
  tf <- tempfile()
  on.exit(unlink(tf))
  utils::download.file(paste0(src, "matrix.bin"), tf, mode = "wb", quiet = TRUE,
                       headers = c(Range = sprintf("bytes=%.0f-%.0f", offset, offset + nbytes - 1)))
  if (file.size(tf) != nbytes)
    stop("the data server ignored the byte range for matrix.bin (got ", file.size(tf),
         " bytes, expected ", nbytes, ")")
  readBin(tf, "raw", nbytes)
}

read_tsv_gz <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  utils::read.delim(con, quote = "", na.strings = "", stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "")
}

open_store <- function(src) {
  if (!is_url(src)) src <- normalizePath(src, mustWork = TRUE)
  else if (!endsWith(src, "/")) src <- paste0(src, "/")
  meta <- as.list(read.dcf(store_file(src, "meta.dcf"))[1, ])
  meta$levels <- strsplit(meta$levels, ",", fixed = TRUE)[[1]]
  cells <- read_tsv_gz(store_file(src, "cells.tsv.gz"))
  genes <- read_tsv_gz(store_file(src, "genes.tsv.gz"))
  list(src = src, meta = meta, cells = cells, genes = genes,
       key = setNames(seq_len(nrow(genes)), toupper(genes$gene)))
}

# Matches user symbols case-insensitively, so "PAX6" and "Pax6" both resolve.
match_genes <- function(store, symbols) {
  hit <- store$key[toupper(symbols)]
  list(found = store$genes$gene[hit[!is.na(hit)]], missing = symbols[is.na(hit)])
}

read_genes <- function(store, genes) {
  rows <- store$key[toupper(genes)]
  blocks <- lapply(rows, function(r)
    decode_gene(store_range(store$src, store$genes$offset[r], store$genes$nbytes[r])))
  setNames(blocks, store$genes$gene[rows])
}

# Per gene x stage x cell type: cell count, % of cells with >= 1 UMI, and mean
# log(1 + 10^4 * count / total) over all cells of the group.
summarise_genes <- function(store, blocks, level, stages) {
  cells <- store$cells
  keep <- cells$stage %in% stages
  grp <- interaction(cells$stage, cells[[level]], drop = TRUE, sep = "\r")
  grp[!keep] <- NA
  grp <- droplevels(grp)
  n <- tabulate(grp, nlevels(grp))
  parts <- do.call(rbind, strsplit(levels(grp), "\r", fixed = TRUE))
  out <- lapply(names(blocks), function(g) {
    b <- blocks[[g]]
    gi <- grp[b$cell]
    ok <- !is.na(gi)
    ln <- log1p(1e4 * b$count[ok] / cells$total[b$cell[ok]])
    data.frame(gene = g, stage = parts[, 1], celltype = parts[, 2], n_cells = n,
               pct_expressing = 100 * tabulate(gi[ok], nlevels(grp)) / n,
               mean_lognorm = vapply(split(ln, factor(as.integer(gi[ok]), seq_len(nlevels(grp)))), sum, 0) / n,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  out$stage_rank <- cells$stage_rank[match(out$stage, cells$stage)]
  out[order(out$stage_rank, out$celltype, match(out$gene, names(blocks))), ]
}

# Per gene pair x stage x cell type: observed co-detection against independence,
# one-sided hypergeometric enrichment test, BH adjustment over the whole table.
coexpression <- function(store, blocks, level, stages) {
  if (length(blocks) < 2) return(NULL)
  cells <- store$cells
  grp <- interaction(cells$stage, cells[[level]], drop = TRUE, sep = "\r")
  grp[!cells$stage %in% stages] <- NA
  grp <- droplevels(grp)
  ng <- nlevels(grp)
  n <- tabulate(grp, ng)
  parts <- do.call(rbind, strsplit(levels(grp), "\r", fixed = TRUE))
  det <- lapply(blocks, function(b) b$cell[!is.na(grp[b$cell])])
  cnt <- lapply(det, function(i) tabulate(grp[i], ng))
  pairs <- utils::combn(names(blocks), 2)
  out <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
    a <- pairs[1, k]; b <- pairs[2, k]
    both <- tabulate(grp[intersect(det[[a]], det[[b]])], ng)
    data.frame(gene1 = a, gene2 = b, stage = parts[, 1], celltype = parts[, 2], n_cells = n,
               n_gene1 = cnt[[a]], n_gene2 = cnt[[b]], n_both = both, stringsAsFactors = FALSE)
  }))
  out$pct_both <- 100 * out$n_both / out$n_cells
  out$expected_both <- out$n_gene1 * out$n_gene2 / out$n_cells
  out$log2_obs_exp <- log2((out$n_both + 0.5) / (out$expected_both + 0.5))
  union <- out$n_gene1 + out$n_gene2 - out$n_both
  out$jaccard <- ifelse(union > 0, out$n_both / union, 0)
  out$p_enrich <- stats::phyper(out$n_both - 1, out$n_gene1, out$n_cells - out$n_gene1,
                                out$n_gene2, lower.tail = FALSE)
  out$p_adj <- stats::p.adjust(out$p_enrich, "BH")
  out
}
