#!/usr/bin/env Rscript
# Builds the "souris_dev" store: mouse forebrain development from E9 to adult,
# merged from eight public studies, with one harmonised column each for age,
# region and cell type (extract/celltypes.R), plus the study, the published
# annotation and whether the cells were sorted on a reporter.
#
# Every study keeps at most CAP cells per stage x region x harmonised type
# (fixed seed), so rare types survive and abundant ones do not dominate.
#
# Usage: Rscript extract/mouse_dev.R <raw_dir> <out_dir>
#   <raw_dir> holds dev_all.loom, l5_all_rev1.loom, geo/ (extract/fetch_mouse_dev.sh)
#   and abc/ (WMB-10Xv3-STR-raw.h5ad, WMB-10Xv3-HPF-raw.h5ad,
#   cell_metadata_with_cluster_annotation.csv).

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "..", "app", "store.R"))
source(file.path(here, "lib.R"))
source(file.path(here, "celltypes.R"))

CAP <- 200L
SEED <- 1L
args <- commandArgs(TRUE)
raw <- args[1]; out <- file.path(args[2], "souris_dev")
geo <- file.path(raw, "geo")

# Days since conception: E<x> = x, P<x> = 19.5 + x (birth at about E19.5).
age_rank <- function(label) {
  x <- as.numeric(sub("^[EP]", "", label))
  ifelse(startsWith(label, "E"), x, 19.5 + x)
}
fmt_e <- function(x) sprintf("E%s", trimws(formatC(x, format = "fg", digits = 3)))

# One part = cells x genes counts plus per-cell columns. `strata` decides the
# sampling groups; the harmonised type is computed from the kept cells' markers.
part <- function(counts, cells, study, sorted = "Toutes les cellules", published = NA_character_) {
  cells$total <- rowSums(counts)
  ok <- cells$total >= 500
  counts <- counts[ok, , drop = FALSE]; cells <- cells[ok, , drop = FALSE]
  cells$type <- classify_cells(as.matrix(counts[, intersect(ALL_MARKERS, colnames(counts)), drop = FALSE]),
                               cells$total)
  keep <- sample_groups(paste(cells$stage, cells$region, cells$type), CAP, SEED)
  cells <- cells[keep, , drop = FALSE]
  cells$study <- study
  cells$sorted <- if (length(sorted) == 1) sorted else sorted[ok][keep]
  if (!"published" %in% names(cells)) cells$published <- published
  message(sprintf("  %-28s %6d cells kept (of %d)", study, nrow(cells), sum(ok)))
  list(counts = counts[keep, , drop = FALSE], cells = cells)
}

parts <- list()

# --- La Manno 2021: forebrain E9-E18 (published Class/Subclass) -----------------
local({
  f <- file.path(raw, "dev_all.loom")
  fid <- H5Fopen(f, "H5F_ACC_RDONLY")
  age <- h5_text(fid, "col_attrs/Age"); tis <- h5_text(fid, "col_attrs/Tissue")
  cls <- h5_text(fid, "col_attrs/Class"); sub <- h5_text(fid, "col_attrs/Subclass")
  H5Fclose(fid)
  region <- c(Forebrain = "Prosencéphale entier", ForebrainDorsal = "Pallium dorsal (futur cortex)",
              ForebrainVentral = "Prosencéphale ventral (MGE + LGE + CGE)",
              ForebrainVentroLateral = "Prosencéphale ventro-latéral")
  keep <- which(tis %in% names(region) & cls != "Bad cells")
  # Pre-sample on the published annotation so the full-gene read stays bounded.
  keep <- keep[sample_groups(paste(age[keep], tis[keep], sub[keep]), CAP, SEED)]
  counts <- read_loom_cells(f, keep)
  stage <- fmt_e(as.numeric(sub("^e", "", age[keep])))
  parts$lamanno <<- part(counts, data.frame(stage = stage, region = region[tis[keep]],
                                            region_detail = tis[keep], published = sub[keep]),
                         "La Manno 2021")
})

# --- Zeisel 2018: P12-P60 forebrain regions (published TaxonomyRank4) -----------
local({
  f <- file.path(raw, "l5_all_rev1.loom")
  fid <- H5Fopen(f, "H5F_ACC_RDONLY")
  age <- h5_text(fid, "col_attrs/Age"); tis <- h5_text(fid, "col_attrs/Tissue")
  r4 <- h5_text(fid, "col_attrs/TaxonomyRank4")
  H5Fclose(fid)
  region <- c(Amygd = "Amygdale", SScortex = "Cortex somatosensoriel", Ctx1 = "Cortex (autres aires)",
              Ctx1.5 = "Cortex (autres aires)", Ctx2 = "Cortex (autres aires)", Ctx3 = "Cortex (autres aires)",
              HC = "Hippocampe", CA1 = "Hippocampe", DentGyr = "Hippocampe",
              StriatDor = "Striatum", StriatVent = "Striatum")
  ages <- lapply(strsplit(sub("^p", "", age), ",\\s*p?"), as.numeric)
  # Pooled samples whose ages span more than a week cannot be placed on the time axis.
  span <- vapply(ages, function(a) diff(range(a)), 0)
  keep <- which(tis %in% names(region) & span <= 7)
  keep <- keep[sample_groups(paste(age[keep], tis[keep], r4[keep]), CAP, SEED)]
  counts <- read_loom_cells(f, keep)
  ages <- ages[keep]
  stage <- vapply(ages, function(a) if (length(a) == 1) sprintf("P%g", a) else sprintf("P%g–P%g", min(a), max(a)), "")
  parts$zeisel <<- part(counts, data.frame(stage = stage, stage_rank = 19.5 + vapply(ages, mean, 0),
                                           region = region[tis[keep]], region_detail = tis[keep],
                                           published = r4[keep]),
                        "Zeisel 2018")
})

# --- Di Bella 2021: somatosensory cortex E10.5-P4 -------------------------------
local({
  # Stages from the GEO sample records (file names round some of them).
  st <- c(GSM5277843 = "E10.5", GSM4635072 = "E11.5", GSM4635073 = "E12.5", GSM4635074 = "E13.5",
          GSM4635075 = "E14.5", GSM4635076 = "E15.5", GSM4635077 = "E16.5", GSM5277844 = "E17.5",
          GSM4635078 = "E18.5", GSM4635079 = "E18.5", GSM4635080 = "P1", GSM4635081 = "P1",
          GSM5277845 = "P4")
  files <- file.path(geo, "GSE153162", list.files(file.path(geo, "GSE153162"), "\\.h5$"))
  stage <- unname(st[sub("_.*", "", basename(files))])
  stopifnot(!anyNA(stage))
  ms <- lapply(files, read_10x_h5)
  cells <- data.frame(stage = rep(stage, vapply(ms, nrow, 0L)))
  parts$dibella <<- part(rbind_genes(ms), transform(cells, region = "Cortex somatosensoriel",
                                                    region_detail = "somatosensory cortex"),
                         "Di Bella 2021")
})

# --- Lee 2022: E12.5 cortex, MGE, LGE, CGE dissected separately ------------------
local({
  s <- c(GSM5090774_CTX = "Pallium dorsal (futur cortex)", GSM5090775_MGE = "MGE",
         GSM5090777_LGE = "LGE", GSM5090776_CGE = "CGE")
  ms <- lapply(names(s), function(k) read_10x_mtx(file.path(geo, "GSE167013", paste0(k, "_"))))
  n <- vapply(ms, nrow, 0L)
  cells <- data.frame(region = rep(unname(s), n), region_detail = rep(sub("^GSM[0-9]+_", "", names(s)), n))
  parts$lee <<- part(rbind_genes(ms), transform(cells, stage = "E12.5"), "Lee 2022")
})

# --- Allaway 2021: Dlx6 lineage, E13 MGE -> P28 cortex ---------------------------
local({
  s <- data.frame(
    key = c("GSM5014302_MGE_Dlx6pos", "GSM5014303_MGE_Dlx6neg", "GSM5014304_P2_Cortex_Dlxpos",
            "GSM5014305_P10_ALM_Dlxpos", "GSM5014306_P10_V1_Dlxpos", "GSM5014307_P28_ALM_Dlxpos",
            "GSM5014308_P28_V1_Dlxpos"),
    stage = c("E13", "E13", "P2", "P10", "P10", "P28", "P28"),
    region = c("MGE", "MGE", "Cortex (autres aires)", "Cortex moteur (ALM)", "Cortex visuel (V1)",
               "Cortex moteur (ALM)", "Cortex visuel (V1)"),
    sorted = c("Triées Dlx6+ (lignée GABAergique)", "Triées Dlx6− (MGE)",
               rep("Triées Dlx+ (interneurones)", 5)))
  ms <- lapply(s$key, function(k) read_10x_mtx(file.path(geo, "GSE164570", paste0(k, "_RNA_"))))
  n <- vapply(ms, nrow, 0L)
  cells <- s[rep(seq_len(nrow(s)), n), c("stage", "region")]
  cells$region_detail <- rep(sub("^GSM[0-9]+_", "", s$key), n)
  parts$allaway <<- part(rbind_genes(ms), cells, "Allaway 2021", sorted = rep(s$sorted, n))
})

# --- Mayer 2018: Lhx6 / Dlx6 lineage, E13.5 MGE/CGE -> P10 ----------------------
local({
  m <- read_dense_counts(file.path(geo, "GSE104156", "GSE104156_digital_expression.csv.gz"), ",")
  lib <- sub("_[ACGT]+$", "", rownames(m))
  info <- data.frame(
    lib = c("MGE_E13.5_Lhx6pos", "CGE_E13.5_Lhx6neg", "Cortex_E18.5", "Subcortex_E18.5",
            "Cortex_P10", "Subcortex_P10"),
    stage = c("E13.5", "E13.5", "E18.5", "E18.5", "P10", "P10"),
    region = c("MGE", "CGE", "Cortex (autres aires)", "Sous-cortex (striatum et voisins)",
               "Cortex (autres aires)", "Sous-cortex (striatum et voisins)"),
    sorted = c("Triées Lhx6+ (MGE)", "Triées Lhx6− (CGE)", rep("Triées Dlx6+ (interneurones)", 4)))
  i <- match(lib, info$lib)
  stopifnot(!anyNA(i))
  cells <- data.frame(stage = info$stage[i], region = info$region[i], region_detail = lib)
  parts$mayer <<- part(m, cells, "Mayer 2018", sorted = info$sorted[i])
})

# --- Hochgerner 2018: dentate gyrus E16.5 -> P132 (published clusters) -----------
local({
  m <- read_dense_counts(file.path(geo, "GSE104323", "GSE104323_10X_expression_data_V2.tab.gz"), "\t")
  md <- read.delim(gzfile(file.path(geo, "GSE104323", "GSE104323_metadata_barcodes_24185cells.txt.gz")),
                   check.names = FALSE, stringsAsFactors = FALSE)
  md <- md[match(rownames(m), md[[1]]), ]
  stopifnot(!anyNA(md[[1]]))
  parts$hochgerner <<- part(m, data.frame(stage = md[["characteristics: age"]], region = "Hippocampe",
                                          region_detail = "Gyrus denté",
                                          published = md[["characteristics: cell cluster"]]),
                            "Hochgerner 2018")
})

# --- Yao 2023 (Allen ABC): adult striatum, amygdala (sAMY), hippocampus ----------
local({
  md <- data.table::fread(file.path(raw, "abc", "cell_metadata_with_cluster_annotation.csv"),
                          select = c("cell_label", "feature_matrix_label", "region_of_interest_acronym",
                                     "class", "subclass"))
  region <- c(sAMY = "Amygdale", STRd = "Striatum", STRv = "Striatum", HIP = "Hippocampe")
  for (fm in c("WMB-10Xv3-STR", "WMB-10Xv3-HPF")) {
    f <- file.path(raw, "abc", paste0(fm, "-raw.h5ad"))
    fid <- H5Fopen(f, "H5F_ACC_RDONLY")
    lab <- h5_text(fid, "obs/cell_label"); genes <- h5_text(fid, "var/gene_symbol")
    H5Fclose(fid)
    m <- md[match(lab, md$cell_label)]
    ok <- which(m$region_of_interest_acronym %in% names(region) & !is.na(m$class))
    rows <- ok[sample_groups(paste(m$region_of_interest_acronym[ok], m$subclass[ok]), CAP, SEED)]
    counts <- read_h5ad_rows(f, "X", rows, genes)
    parts[[fm]] <<- part(counts, data.frame(stage = "P56–P69", stage_rank = 19.5 + 60,
                                            region = region[m$region_of_interest_acronym[rows]],
                                            region_detail = m$region_of_interest_acronym[rows],
                                            published = m$subclass[rows]),
                         "Yao 2023 (Allen ABC)")
  }
})

# --- merge and write ---------------------------------------------------------------
cells <- do.call(rbind, lapply(parts, function(p) {
  c <- p$cells
  if (!"stage_rank" %in% names(c)) c$stage_rank <- age_rank(c$stage)
  c[, c("stage", "stage_rank", "region", "region_detail", "type", "study", "sorted", "published")]
}))
counts <- rbind_genes(lapply(parts, `[[`, "counts"))

# Time windows around the requested ages (days since conception, birth = 19.5),
# cut half-way between neighbouring targets. P10 has data of its own and keeps
# its own window. meta$stage_note lists the real ages inside each window.
targets <- c("≈ E10.5" = 10.5, "≈ E12.5" = 12.5, "≈ E13.5" = 13.5, "≈ E15.5" = 15.5, "≈ P0" = 19.5,
             "≈ P10" = 29.5, "≈ P21" = 40.5, "≈ P40" = 59.5, "≈ P100" = 119.5)
cuts <- c(-Inf, head(targets, -1) + diff(targets) / 2)
w <- findInterval(cells$stage_rank, cuts)
real <- tapply(seq_along(w), w, function(i) {
  s <- unique(cells$stage[i][order(cells$stage_rank[i])])
  paste(s, collapse = ", ")
})
cells$age_exact <- cells$stage
cells$stage <- names(targets)[w]
cells$stage_rank <- unname(targets[w])
stage_note <- paste(sprintf("%s : %s", names(targets)[as.integer(names(real))], real), collapse = " ; ")

ord <- order(cells$stage_rank, cells$region, cells$study)
cells <- cells[ord, ]; counts <- counts[ord, ]
names(cells)[match(c("type", "region", "study", "published", "age_exact", "sorted"), names(cells))] <-
  c("Type", "Region", "Study", "Published", "Age", "Sorting")
write_store_matrix(out, cells = cells, counts = counts, meta = list(
  id = "souris_dev", species = "mouse",
  label = "Souris : développement du prosencéphale, E9 → adulte (8 études)",
  citation = paste("La Manno 2021; Zeisel 2018; Di Bella 2021; Lee 2022; Allaway 2021; Mayer 2018;",
                   "Hochgerner 2018; Yao 2023 (Allen Brain Cell Atlas)"),
  levels = c("Type", "Region", "Study", "Age", "Sorting"),
  filter_level = "Region",
  stage_note = stage_note,
  sampling = sprintf("au plus %d cellules par étude × stade × région × type, tirage aléatoire (graine %d)",
                     CAP, SEED)))
