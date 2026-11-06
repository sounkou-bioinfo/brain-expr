# ============================================================================
# CONFIGURATION : la seule partie à modifier
# ============================================================================

# Les ~10 gènes d'intérêt (symboles souris, ex. "Pax6", pas "PAX6").
# Liste d'EXEMPLE (marqueurs classiques) : à remplacer par la sienne.
genes <- c("Pax6", "Sox2", "Hes5", "Mki67", "Eomes",
           "Neurod2", "Tbr1", "Gad2", "Olig2", "Aqp4")

# Effectif minimal d'un groupe (stade x type cellulaire) pour être rapporté
min_cells <- 20

datasets <- list(
  # La Manno et al. 2021 (Nature) : cerveau embryonnaire, E7 -> E18, un jour sur un
  lamanno2021 = list(
    file   = "data/raw/dev_all.loom",          # ~19 Go, mousebrain.org/development
    stages = c("E12.5", "E13.5"),
    stage_window = 0,                           # jours ; 0 = stade exact
    age_attrs = c("Age"),
    numeric_age_prefix = "E",                   # si Age est stocké en 12.5 sans préfixe
    # niveaux d'annotation, du plus grossier au plus fin ; ceux absents sont ignorés
    celltype_attrs = c("Class", "Subclass", "ClusterName"),
    tissue_attrs = c("Tissue", "Region"),
    tissue_exclude = NULL
  ),
  # Zeisel et al. 2018 (Cell) : système nerveux adolescent, P12 -> P30 (+ 6 et 8 sem.)
  zeisel2018 = list(
    file   = "data/raw/l5_all.loom",           # mousebrain.org/adolescent
    stages = c("P21"),
    stage_window = 2,                           # P19-P23 : peu de cellules à P21 pile
    age_attrs = c("Age"),
    numeric_age_prefix = "P",
    celltype_attrs = c("TaxonomyRank1", "TaxonomyRank3", "TaxonomyRank4",
                       "Class", "Subclass", "ClusterName"),
    tissue_attrs = c("Tissue", "Region"),
    # cerveau seulement : retirer ganglions, entérique, sympathique, moelle
    tissue_exclude = "^(DRG|ENS|Enteric|Sympath|SpinalC|SpinalCord|Spinal)"
  )
)
