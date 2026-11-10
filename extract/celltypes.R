# Harmonised cell types for the mouse atlases, assigned from marker genes with
# one rule for every study, so the same label means the same thing across ages
# and datasets. Published annotations are kept in their own column.
#
# Each class has a marker set; a cell's class score is the mean of
# log(1 + 10^4 * count / total) over the markers. The best-scoring class wins
# when its score reaches MIN_SCORE, otherwise the cell is "Non classé".
# GABAergic neurons and progenitors are then split on single markers.

MIN_SCORE <- 0.5

CLASS_MARKERS <- list(
  progenitor = c("Mki67", "Top2a", "Nes", "Hes5", "Fabp7", "Sox2", "Ascl1", "Eomes", "Neurog2"),
  gaba = c("Gad1", "Gad2", "Slc32a1", "Dlx1", "Dlx2", "Dlx5", "Dlx6os1"),
  glut = c("Neurod2", "Neurod6", "Neurod1", "Slc17a7", "Slc17a6", "Tbr1", "Satb2", "Prox1", "Trp73"),
  astro = c("Aqp4", "Aldh1l1", "Gja1", "Agt", "Slc1a2"),
  oligo = c("Olig1", "Sox10", "Pdgfra", "Mbp", "Plp1", "Mog"),
  immune = c("C1qa", "C1qb", "Cx3cr1", "P2ry12", "Csf1r", "Tyrobp", "Mrc1"),
  vascular = c("Cldn5", "Flt1", "Pecam1", "Pdgfrb", "Rgs5", "Col1a2", "Dcn", "Hba-a1", "Hbb-bs")
)
SUBTYPE_MARKERS <- c("Nkx2-1", "Lhx6", "Pvalb", "Sst", "Ppp1r1b", "Penk", "Tac1")

TYPE_LEVELS <- c(
  "Progéniteurs Nkx2-1+ (MGE)", "Autres progéniteurs (glie radiaire, IPC)",
  "Neurones GABAergiques Dlx+ (post-mitotiques, autres)",
  "Interneurones Lhx6+ immatures (issus de la MGE)",
  "Interneurones PV+", "Interneurones SST+",
  "Neurones de projection du striatum (MSN)",
  "Neurones excitateurs (pyramidaux, granulaires)",
  "Astrocytes", "Lignée oligodendrocytaire (OPC, oligodendrocytes)",
  "Microglie / macrophages", "Cellules vasculaires, méninges, sang", "Non classé")

ALL_MARKERS <- unique(c(unlist(CLASS_MARKERS), SUBTYPE_MARKERS))

# counts: dense cells x genes matrix of raw counts with column names (missing
# markers are ignored); total: total UMI per cell.
classify_cells <- function(counts, total) {
  ln <- log1p(1e4 * counts / pmax(total, 1))
  score <- vapply(CLASS_MARKERS, function(m) {
    m <- intersect(m, colnames(ln))
    if (!length(m)) return(rep(0, nrow(ln)))
    rowMeans(ln[, m, drop = FALSE])
  }, numeric(nrow(ln)))
  if (is.null(dim(score))) score <- matrix(score, nrow = 1, dimnames = list(NULL, names(CLASS_MARKERS)))
  best <- colnames(score)[max.col(score, ties.method = "first")]
  best[apply(score, 1, max) < MIN_SCORE] <- "none"
  det <- function(g) if (g %in% colnames(counts)) counts[, g] > 0 else rep(FALSE, nrow(counts))
  msn <- if (all(c("Ppp1r1b", "Penk", "Tac1") %in% colnames(ln)))
    rowMeans(ln[, c("Ppp1r1b", "Penk", "Tac1")]) >= 1 else rep(FALSE, nrow(ln))
  out <- rep("Non classé", nrow(counts))
  out[best == "progenitor"] <- ifelse(det("Nkx2-1")[best == "progenitor"],
                                      TYPE_LEVELS[1], TYPE_LEVELS[2])
  g <- best == "gaba"
  out[g] <- TYPE_LEVELS[3]
  out[g & det("Lhx6")] <- TYPE_LEVELS[4]
  out[g & msn & !det("Lhx6")] <- TYPE_LEVELS[7]
  out[g & det("Sst")] <- TYPE_LEVELS[6]
  out[g & det("Pvalb")] <- TYPE_LEVELS[5]
  out[best == "glut"] <- TYPE_LEVELS[8]
  out[best == "astro"] <- TYPE_LEVELS[9]
  out[best == "oligo"] <- TYPE_LEVELS[10]
  out[best == "immune"] <- TYPE_LEVELS[11]
  out[best == "vascular"] <- TYPE_LEVELS[12]
  out
}
