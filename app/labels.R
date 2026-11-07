# Readable French names for the annotation levels, codes and stages of each atlas.
# Cell type names that are standard English terms (Radial glia, Neuroblast...)
# stay as published; cryptic codes and abbreviations are spelled out.

LEVEL_LABELS <- list(
  lamanno2021 = c(Class = "Classe cellulaire", Subclass = "Sous-classe (plus fin)",
                  Tissue = "Région disséquée"),
  zeisel2018 = c(TaxonomyRank3 = "Grande famille", TaxonomyRank4 = "Type cellulaire (plus fin)",
                 Tissue = "Région disséquée"),
  braun2023 = c(CellClass = "Classe cellulaire", Region = "Région disséquée"),
  velmeshev2023 = c(Lineage = "Lignée cellulaire", Region = "Région du cortex",
                    Study = "Étude d'origine")
)

VALUE_LABELS <- list(
  lamanno2021 = list(Tissue = c(
    ForebrainDorsal = "Prosencéphale dorsal", ForebrainVentral = "Prosencéphale ventral",
    Midbrain = "Mésencéphale", Hindbrain = "Rhombencéphale")),
  zeisel2018 = list(Tissue = c(
    Amygd = "Amygdale", CA1 = "Hippocampe (CA1)", HC = "Hippocampe", CB = "Cervelet",
    Ctx1 = "Cortex (lot 1)", Ctx2 = "Cortex (lot 2)", Ctx3 = "Cortex (lot 3)",
    Hypoth = "Hypothalamus", Thal = "Thalamus", MBd = "Mésencéphale dorsal",
    MBv = "Mésencéphale ventral", Pons = "Pont", Medulla = "Bulbe rachidien")),
  braun2023 = list(Region = c(
    Brain = "Cerveau entier", Head = "Tête entière", Forebrain = "Prosencéphale",
    Telencephalon = "Télencéphale", Diencephalon = "Diencéphale", Midbrain = "Mésencéphale",
    Hindbrain = "Rhombencéphale", Cerebellum = "Cervelet", Pons = "Pont",
    Medulla = "Bulbe rachidien")),
  velmeshev2023 = list(
    Lineage = c(ExNeu = "Neurones excitateurs", IN = "Interneurones (inhibiteurs)",
                AST = "Astrocytes", OPC = "Précurseurs d'oligodendrocytes (OPC)",
                OL = "Oligodendrocytes", MG = "Microglie", VASC = "Cellules vasculaires",
                GLIALPROG = "Progéniteurs gliaux", OUT = "Autres / non classées"),
    Region = c(FC = "Cortex frontal", CC = "Cortex cingulaire", TC = "Cortex temporal",
               IC = "Cortex insulaire", MC = "Cortex moteur", CTX = "Cortex (non précisé)",
               GE = "Éminence ganglionnaire", MGE = "Éminence ganglionnaire médiane",
               LGE = "Éminence ganglionnaire latérale", CGE = "Éminence ganglionnaire caudale"))
)

STAGE_LABELS <- list(
  lamanno2021 = function(x) x,
  zeisel2018 = function(x) "P19–P23 (≈ 3 semaines après la naissance)",
  braun2023 = function(x) paste(sub(".", ",", sub(" pcw$", "", x), fixed = TRUE), "sem."),
  velmeshev2023 = function(x) unname(c(
    "2nd trimester" = "2e trimestre", "3rd trimester" = "3e trimestre",
    "0-1 years" = "0–1 an", "1-2 years" = "1–2 ans", "2-4 years" = "2–4 ans",
    "4-10 years" = "4–10 ans", "10-20 years" = "10–20 ans", "Adult" = "Adulte")[x])
)

# Plain description of each atlas for the guide and the header.
ATLAS_TEXT <- list(
  lamanno2021 = list(
    short = "Souris · embryon E12.5–E13.5",
    what = "Cerveau de souris embryonnaire, disséqué par région, aux jours E12.5 et E13.5 : la période où la neurogenèse bat son plein.",
    stages = "E12.5 et E13.5 = jours de gestation."),
  zeisel2018 = list(
    short = "Souris · P19–P23",
    what = "Système nerveux de souris juvéniles (P19 à P23), ici réduit au cerveau : moelle et ganglions sont retirés.",
    stages = "Un seul stade : P19–P23 = 19 à 23 jours après la naissance."),
  braun2023 = list(
    short = "Humain · 5–14 semaines",
    what = "Cerveau humain du premier trimestre, de 5 à 14 semaines après la conception, disséqué par région.",
    stages = "« 6,6 sem. » = 6,6 semaines après la conception (pcw en anglais)."),
  velmeshev2023 = list(
    short = "Humain · 2e trimestre → adulte",
    what = "Cortex humain du deuxième trimestre de grossesse jusqu'à l'âge adulte ; réunit quatre études (Velmeshev, Ramos, Herring, Trevino).",
    stages = "Tranches d'âge, du 2e trimestre de grossesse à l'adulte.")
)

# Replaces codes and stages in a store's cell table by their readable names.
relabel_store <- function(s) {
  id <- s$meta$id
  for (lv in names(VALUE_LABELS[[id]])) {
    v <- s$cells[[lv]]
    hit <- v %in% names(VALUE_LABELS[[id]][[lv]])
    v[hit] <- VALUE_LABELS[[id]][[lv]][v[hit]]
    s$cells[[lv]] <- v
  }
  s$cells$stage <- STAGE_LABELS[[id]](s$cells$stage)
  s$level_names <- setNames(s$meta$levels, LEVEL_LABELS[[id]][s$meta$levels])
  s$text <- ATLAS_TEXT[[id]]
  s
}
