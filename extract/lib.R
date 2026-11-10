# Shared extraction helpers: store writer and readers for the input formats
# (loom, AnnData h5ad, 10x mtx and h5, dense count tables).
#
# Readers return cells x genes count matrices (Matrix::dgCMatrix) whose column
# names are gene symbols; duplicated symbols are summed.

suppressPackageStartupMessages({ library(rhdf5); library(Matrix) })

h5_text <- function(fid, path) as.character(h5read(fid, path))

# AnnData categorical column: integer codes into a categories array.
h5ad_factor <- function(fid, path) {
  x <- h5read(fid, path)
  if (!is.list(x)) return(as.character(x))
  as.character(x$categories)[x$codes + 1L]
}

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
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

# Writes a store from a cells x genes dgCMatrix: each column is one gene block.
write_store_matrix <- function(out, meta, cells, counts) {
  cells$total <- rowSums(counts)
  write_store(out, meta, cells, colnames(counts), function(emit) {
    p <- counts@p
    for (g in seq_len(ncol(counts))) {
      k <- if (p[g + 1] > p[g]) (p[g] + 1):p[g + 1] else integer()
      emit(g, counts@i[k], as.integer(round(counts@x[k])))
    }
  })
}

# Sums columns that share a gene symbol.
collapse_symbols <- function(m, symbols) {
  if (!anyDuplicated(symbols)) { colnames(m) <- symbols; return(m) }
  f <- factor(symbols, unique(symbols))
  m <- m %*% Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                                  dims = c(length(f), nlevels(f)))
  colnames(m) <- levels(f)
  as(m, "CsparseMatrix")
}

# Binds cells x genes matrices with different gene sets on the union of genes.
rbind_genes <- function(ms) {
  genes <- unique(unlist(lapply(ms, colnames)))
  ms <- lapply(ms, function(m) {
    t <- as(m, "TsparseMatrix")
    Matrix::sparseMatrix(i = t@i + 1L, j = match(colnames(m), genes)[t@j + 1L], x = t@x,
                         dims = c(nrow(m), length(genes)), dimnames = list(NULL, genes))
  })
  do.call(rbind, ms)
}

# 10x Matrix Market triplet: <prefix>matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz.
read_10x_mtx <- function(prefix) {
  m <- Matrix::readMM(gzfile(paste0(prefix, "matrix.mtx.gz")))
  f <- read.delim(gzfile(paste0(prefix, "features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  b <- readLines(gzfile(paste0(prefix, "barcodes.tsv.gz")))
  m <- Matrix::t(as(m, "CsparseMatrix"))
  rownames(m) <- b
  collapse_symbols(as(m, "CsparseMatrix"), f[[2]])
}

# 10x HDF5 (Cell Ranger v2 "genome/gene_names" or v3 "matrix/features/name").
read_10x_h5 <- function(file) {
  fid <- H5Fopen(file, "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))
  grp <- h5ls(fid, recursive = FALSE)$name[1]
  sym <- if (grp == "matrix") h5_text(fid, "matrix/features/name") else h5_text(fid, paste0(grp, "/gene_names"))
  shape <- as.integer(h5read(fid, paste0(grp, "/shape")))
  m <- Matrix::sparseMatrix(i = as.integer(h5read(fid, paste0(grp, "/indices"))) + 1L,
                            p = as.integer(h5read(fid, paste0(grp, "/indptr"))),
                            x = as.numeric(h5read(fid, paste0(grp, "/data"))), dims = shape)
  m <- Matrix::t(m)
  rownames(m) <- h5_text(fid, paste0(grp, "/barcodes"))
  collapse_symbols(as(m, "CsparseMatrix"), sym)
}

# Dense genes x cells table (first column = gene symbol), read in row chunks.
read_dense_counts <- function(file, sep, chunk = 2000L) {
  con <- gzfile(file, "rt")
  on.exit(close(con))
  header <- scan(text = readLines(con, n = 1), what = "", sep = sep, quote = "\"", quiet = TRUE)
  cells <- header[-1]
  parts <- list(); genes <- character()
  repeat {
    lines <- readLines(con, n = chunk)
    if (!length(lines)) break
    tab <- data.table::fread(text = lines, sep = sep, header = FALSE, quote = "\"")
    genes <- c(genes, tab[[1]])
    parts[[length(parts) + 1L]] <- as(as.matrix(tab[, -1]), "CsparseMatrix")
  }
  m <- Matrix::t(do.call(rbind, parts))
  rownames(m) <- cells
  collapse_symbols(as(m, "CsparseMatrix"), genes)
}

# Selected rows of an AnnData CSR matrix (cells x genes), as a dgCMatrix.
read_h5ad_rows <- function(file, xgrp, rows, genes, batch = 2000L) {
  fid <- H5Fopen(file, "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))
  indptr <- h5read(fid, paste0(xgrp, "/indptr"), bit64conversion = "double")
  d <- H5Dopen(fid, paste0(xgrp, "/data")); di <- H5Dopen(fid, paste0(xgrp, "/indices"))
  on.exit({ H5Dclose(d); H5Dclose(di) }, add = TRUE, after = FALSE)
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
    out <- Matrix::sparseMatrix(i = rep.int(seq_along(k), n),
                                j = as.integer(H5Dread(di, h5spaceFile = spi, h5spaceMem = ms)) + 1L,
                                x = as.numeric(H5Dread(d, h5spaceFile = sp, h5spaceMem = ms)),
                                dims = c(length(k), length(genes)))
    H5Sclose(ms); H5Sclose(sp); H5Sclose(spi)
    out
  })
  collapse_symbols(do.call(rbind, parts), genes)
}

# Loom (genes x cells) restricted to `keep` cells, as a cells x genes dgCMatrix.
read_loom_cells <- function(file, keep, band = 512L) {
  fid <- H5Fopen(file, "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))
  genes <- h5_text(fid, "row_attrs/Gene")
  parts <- lapply(seq(1L, length(genes), by = band), function(s) {
    g <- s:min(length(genes), s + band - 1L)
    as(h5read(fid, "matrix", index = list(NULL, g))[keep, , drop = FALSE], "CsparseMatrix")
  })
  collapse_symbols(do.call(cbind, parts), genes)
}

# Fixed-seed sample of at most `cap` cells per group; returns sorted indices.
sample_groups <- function(group, cap, seed = 1L) {
  set.seed(seed)
  idx <- split(seq_along(group), group, drop = TRUE)
  sort(unlist(lapply(idx, function(i) if (length(i) > cap) i[sample.int(length(i), cap)] else i),
              use.names = FALSE))
}
