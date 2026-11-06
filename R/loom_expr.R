# Fonctions pour extraire l'expression de quelques gènes depuis des fichiers
# .loom (mousebrain.org) sans jamais charger la matrice entière en mémoire.
# Dépendances : hdf5r, data.table

suppressPackageStartupMessages({
  library(hdf5r)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Âges : "e12.5", "E12.5", "12.5", "p21", "P21", "6w" -> stade normalisé
# ---------------------------------------------------------------------------
parse_age <- function(x, numeric_prefix = "E") {
  x <- trimws(tolower(as.character(x)))
  out <- data.table(raw = x, prefix = NA_character_, day = NA_real_)
  m <- regmatches(x, regexec("^([ep])\\s*([0-9]+(\\.[0-9]+)?)$", x))
  ok <- lengths(m) > 0
  out[ok, `:=`(prefix = toupper(vapply(m[ok], `[`, "", 2)),
               day = as.numeric(vapply(m[ok], `[`, "", 3)))]
  wk <- regmatches(x, regexec("^([0-9]+)\\s*w(k|eeks?)?$", x))
  okw <- lengths(wk) > 0
  out[okw, `:=`(prefix = "P", day = 7 * as.numeric(vapply(wk[okw], `[`, "", 2)))]
  num <- !ok & !okw & grepl("^[0-9]+(\\.[0-9]+)?$", x)
  out[num, `:=`(prefix = toupper(numeric_prefix), day = as.numeric(x[num]))]
  out[, label := ifelse(is.na(day), raw,
                        paste0(prefix, as.character(day)))]
  out[]
}

# stage "E12.5" -> list(prefix="E", day=12.5)
parse_stage <- function(s) {
  p <- parse_age(s)
  list(prefix = p$prefix, day = p$day)
}

# ---------------------------------------------------------------------------
# Accès loom
# ---------------------------------------------------------------------------
open_loom <- function(path) {
  if (!file.exists(path)) stop("Fichier introuvable : ", path,
                               "\n  -> téléchargez-le depuis mousebrain.org (voir README)")
  H5File$new(path, mode = "r")
}

loom_attr_names <- function(f, which = c("col_attrs", "row_attrs")) {
  which <- match.arg(which)
  names(f[[which]])
}

read_attr <- function(f, which, name) f[[paste0(which, "/", name)]][]

pick_attr <- function(f, candidates, which = "col_attrs", required = TRUE) {
  have <- loom_attr_names(f, which)
  hit <- candidates[candidates %in% have]
  if (length(hit)) return(hit[1])
  if (required) stop("Aucun des attributs ", paste(candidates, collapse = ", "),
                     " dans ", which, ". Disponibles : ", paste(have, collapse = ", "))
  NULL
}

# Orientation : hdf5r présente la matrice dans l'ordre inverse de h5py/loompy,
# donc en général cellules x gènes. On le détecte au lieu de le supposer.
matrix_layout <- function(f, n_genes, n_cells) {
  d <- f[["matrix"]]$dims
  if (d[1] == n_cells && d[2] == n_genes) return("cells_x_genes")
  if (d[1] == n_genes && d[2] == n_cells) return("genes_x_cells")
  stop("Dimensions de /matrix (", paste(d, collapse = " x "),
       ") incompatibles avec ", n_genes, " gènes et ", n_cells, " cellules")
}

# Lit les comptes d'un gène (toutes les cellules) -> vecteur
read_gene <- function(f, layout, gi) {
  m <- f[["matrix"]]
  if (layout == "cells_x_genes") as.numeric(m[, gi]) else as.numeric(m[gi, ])
}

# Total de comptes par cellule pour un sous-ensemble de cellules, par blocs
# (utilisé seulement si le loom n'a pas d'attribut de total).
compute_totals <- function(f, layout, cells, block = 2000L, verbose = TRUE) {
  m <- f[["matrix"]]
  cells <- sort(unique(cells))
  tot <- numeric(length(cells))
  starts <- seq(1L, length(cells), by = block)
  for (k in seq_along(starts)) {
    idx <- cells[starts[k]:min(starts[k] + block - 1L, length(cells))]
    span <- min(idx):max(idx)
    x <- if (layout == "cells_x_genes") m[span, , drop = FALSE] else t(m[, span, drop = FALSE])
    s <- rowSums(x)[idx - min(idx) + 1L]
    tot[starts[k] - 1L + seq_along(idx)] <- s
    if (verbose && k %% 10 == 0) message(sprintf("  totaux : bloc %d/%d", k, length(starts)))
  }
  setNames(tot, cells)
}

# ---------------------------------------------------------------------------
# Extraction d'un jeu de données
# ---------------------------------------------------------------------------
#' @param cfg liste : file, stages, stage_window (jours), celltype_attrs,
#'   age_attrs, tissue_attrs, tissue_exclude (regex), numeric_age_prefix
#' @return list(cells = data.table métadonnées, expr = data.table long cellule×gène)
extract_dataset <- function(name, cfg, genes, verbose = TRUE) {
  if (verbose) message("== ", name, " : ", cfg$file)
  f <- open_loom(cfg$file)
  on.exit(f$close_all(), add = TRUE)

  gene_attr <- pick_attr(f, c("Gene", "gene", "Genes", "var_names"), "row_attrs")
  all_genes <- as.character(read_attr(f, "row_attrs", gene_attr))

  age_attr <- pick_attr(f, cfg$age_attrs %||% c("Age", "age", "Timepoint"))
  ages <- parse_age(read_attr(f, "col_attrs", age_attr),
                    numeric_prefix = cfg$numeric_age_prefix %||% "E")
  n_cells <- nrow(ages)
  layout <- matrix_layout(f, length(all_genes), n_cells)

  # --- sélection des cellules par stade (avec fenêtre en jours) ---
  win <- cfg$stage_window %||% 0
  stage_of <- rep(NA_character_, n_cells)
  for (s in cfg$stages) {
    st <- parse_stage(s)
    hit <- which(ages$prefix == st$prefix & abs(ages$day - st$day) <= win + 1e-9)
    stage_of[hit] <- toupper(s)
  }
  keep <- !is.na(stage_of)

  # --- filtre tissus (ex. retirer SNP/ENS/moelle pour Zeisel) ---
  tissue_attr <- pick_attr(f, cfg$tissue_attrs %||% c("Tissue", "Region"), required = FALSE)
  tissue <- if (!is.null(tissue_attr)) as.character(read_attr(f, "col_attrs", tissue_attr)) else rep(NA, n_cells)
  if (!is.null(cfg$tissue_exclude) && !is.null(tissue_attr))
    keep <- keep & !grepl(cfg$tissue_exclude, tissue)
  cells <- which(keep)
  if (!length(cells)) {
    warning(name, " : aucune cellule pour les stades ", paste(cfg$stages, collapse = ", "),
            ". Âges présents : ", paste(head(unique(ages$label), 30), collapse = ", "))
    return(NULL)
  }

  # --- annotations de types cellulaires (plusieurs niveaux possibles) ---
  ct_attrs <- intersect(cfg$celltype_attrs, loom_attr_names(f, "col_attrs"))
  if (!length(ct_attrs)) stop(name, " : aucun attribut de type cellulaire parmi ",
                              paste(cfg$celltype_attrs, collapse = ", "))
  meta <- data.table(cell = cells, dataset = name, stage = stage_of[cells],
                     age = ages$label[cells], tissue = tissue[cells])
  for (a in ct_attrs) set(meta, j = paste0("ct_", a),
                          value = as.character(read_attr(f, "col_attrs", a))[cells])

  # --- totaux par cellule pour la normalisation ---
  tot_attr <- pick_attr(f, c("TotalUMI", "TotalMolecules", "nUMI", "total_counts", "n_counts"),
                        required = FALSE)
  if (!is.null(tot_attr)) {
    meta[, total := as.numeric(read_attr(f, "col_attrs", tot_attr))[cells]]
  } else {
    if (verbose) message("  pas d'attribut de total : calcul sur ", length(cells), " cellules")
    tt <- compute_totals(f, layout, cells, verbose = verbose)
    meta[, total := tt[as.character(cell)]]
  }

  # --- gènes demandés ---
  gi <- match(genes, all_genes)
  if (anyNA(gi)) {
    miss <- genes[is.na(gi)]
    # tentative insensible à la casse (Pax6 / PAX6)
    ci <- match(tolower(miss), tolower(all_genes))
    gi[is.na(gi)] <- ci
    still <- genes[is.na(gi)]
    if (length(still)) warning(name, " : gènes absents : ", paste(still, collapse = ", "))
  }
  expr <- rbindlist(lapply(which(!is.na(gi)), function(k) {
    v <- read_gene(f, layout, gi[k])[cells]
    data.table(cell = cells, gene = genes[k], counts = v)
  }))
  expr[meta, total := i.total, on = "cell"]
  expr[, lognorm := log1p(1e4 * counts / pmax(total, 1))]
  expr[, total := NULL]
  if (verbose) message(sprintf("  %d cellules, %d gènes, stades : %s", length(cells),
                               uniqueN(expr$gene), paste(unique(meta$stage), collapse = ", ")))
  list(cells = meta, expr = expr)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Résumés
# ---------------------------------------------------------------------------
summarise_expression <- function(ds, level, min_cells = 20) {
  col <- paste0("ct_", level)
  if (!col %in% names(ds$cells)) return(NULL)
  m <- ds$cells[, .(cell, dataset, stage, celltype = get(col))]
  x <- ds$expr[m, on = "cell"]
  s <- x[, .(n_cells = .N,
             mean_lognorm = mean(lognorm),
             mean_lognorm_expressing = if (any(counts > 0)) mean(lognorm[counts > 0]) else 0,
             pct_expressing = 100 * mean(counts > 0),
             mean_counts = mean(counts)),
         by = .(dataset, stage, celltype, gene)]
  s <- s[n_cells >= min_cells]
  s[, level := level]
  # z-score par gène au sein de chaque jeu de données (comparaisons relatives)
  s[, z_within_dataset := {
    sd_ <- sd(mean_lognorm); if (is.na(sd_) || sd_ == 0) 0 else (mean_lognorm - mean(mean_lognorm)) / sd_
  }, by = .(dataset, gene)]
  setcolorder(s, c("dataset", "level", "stage", "celltype", "gene"))
  s[]
}

# Coexpression par paire de gènes dans chaque groupe (stade x type cellulaire).
# observé = fraction de cellules où les deux gènes sont détectés ;
# attendu = pA * pB (indépendance) ; test hypergéométrique unilatéral.
summarise_coexpression <- function(ds, level, min_cells = 20) {
  col <- paste0("ct_", level)
  if (!col %in% names(ds$cells)) return(NULL)
  m <- ds$cells[, .(cell, dataset, stage, celltype = get(col))]
  w <- dcast(ds$expr[, .(cell, gene, det = counts > 0)], cell ~ gene, value.var = "det")
  w <- w[m, on = "cell"]
  genes <- setdiff(names(w), c("cell", "dataset", "stage", "celltype"))
  if (length(genes) < 2) return(NULL)
  pairs <- t(combn(sort(genes), 2))
  res <- w[, {
    n <- .N
    if (n < min_cells) NULL else {
      M <- as.matrix(.SD); storage.mode(M) <- "integer"
      both <- crossprod(M)              # n_both pour toutes les paires
      nA <- diag(both)
      out <- data.table(gene1 = pairs[, 1], gene2 = pairs[, 2])
      out[, `:=`(n_cells = n,
                 n_gene1 = nA[gene1], n_gene2 = nA[gene2],
                 n_both = both[cbind(gene1, gene2)])]
      out
    }
  }, by = .(dataset, stage, celltype), .SDcols = genes]
  if (!nrow(res)) return(NULL)
  res[, `:=`(pct_both = 100 * n_both / n_cells,
             expected_both = n_gene1 * n_gene2 / n_cells)]
  res[, log2_obs_exp := log2((n_both + 0.5) / (expected_both + 0.5))]
  res[, jaccard := fifelse(n_gene1 + n_gene2 - n_both > 0,
                           n_both / (n_gene1 + n_gene2 - n_both), 0)]
  res[, p_enrich := phyper(n_both - 1, n_gene1, n_cells - n_gene1, n_gene2, lower.tail = FALSE)]
  res[, p_adj := p.adjust(p_enrich, "BH"), by = dataset]
  res[, level := level]
  setcolorder(res, c("dataset", "level", "stage", "celltype", "gene1", "gene2"))
  res[]
}
