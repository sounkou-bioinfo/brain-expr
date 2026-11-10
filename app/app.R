# Explorateur d'expression : n'importe quelle liste de gènes, par type cellulaire et par
# stade, dans quatre atlas scRNA-seq du cerveau en développement (souris et humain).
#
# Les données sont des « stores » par atlas (format : store.R), lus gène par gène.
# - R natif : dossier local, BRAIN_EXPR_DATA ou ../stores à côté de app/.
# - shinylive / app de bureau : dossier data/ servi à côté de l'app, lu par requêtes HTTP Range.

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(ggplot2); library(DT)
})
source("store.R", local = TRUE)
source("labels.R", local = TRUE)

ATLASES <- c("souris_dev", "lamanno2021", "zeisel2018", "braun2023", "velmeshev2023")
EXAMPLES <- list(
  rho_core = "Myo9b Trio Rhoa Rac1 Cdc42",
  markers = "Pax6 Sox2 Hes5 Mki67 Eomes Neurod2 Tbr1 Gad2 Olig2 Aqp4",
  epilepsy = "Scn1a Scn8a Kcnq2 Stxbp1 Cdkl5 Gabra1 Depdc5 Tsc2 Arx Foxg1",
  lineage = "Nkx2-1 Lhx6 Dlx1 Dlx2 Dlx5 Cxcr4 Cxcr7 Erbb4 Sst Pvalb Gad1 Mki67"
)
# Rho GTPase regulators from the HGNC gene groups, with MGI mouse orthologues.
RHO <- utils::read.delim("rho_regulators.tsv", stringsAsFactors = FALSE)
PRESETS <- c("Choisir une liste…" = "",
             "Myo9b, Trio, Rhoa, Rac1, Cdc42" = "rho_core",
             "Régulateurs des Rho GTPases : tous" = "rho:all",
             setNames(paste0("rho:", unique(RHO$family)),
                      paste("Régulateurs Rho :", unique(RHO$family))),
             "Lignée des interneurones (MGE → cortex)" = "lineage",
             "Marqueurs classiques" = "markers",
             "Gènes d'épilepsie" = "epilepsy")
MAX_GENES <- 200L

# Gene symbols of a preset list, in the species' nomenclature.
preset_genes <- function(p, species) {
  if (!startsWith(p, "rho:")) return(strsplit(EXAMPLES[[p]], " ")[[1]])
  fam <- sub("^rho:", "", p)
  r <- if (fam == "all") RHO else RHO[RHO$family == fam, ]
  if (identical(species, "human")) r$human else ifelse(nzchar(r$mouse), r$mouse, r$human)
}
# The dotplot and the coexpression stay readable up to this many genes; the
# heatmap shows them all.
MAX_DOT <- 30L
# Groups smaller than this are hidden: with 2 or 3 cells, one detection is 33-50 %.
MIN_CELLS <- 20L
POOLED <- "Stades sélectionnés, regroupés"
VIDEO_URL <- "https://sounkou-bioinfo.github.io/brain-expr/walkthrough.mp4"

local_root <- function() {
  for (d in c(Sys.getenv("BRAIN_EXPR_DATA"), "../stores", "data"))
    if (nzchar(d) && file.exists(file.path(d, ATLASES[1], "meta.dcf"))) return(normalizePath(d))
  NULL
}

# Directory of the page that hosts the app. Under shinylive the app runs at a
# virtual path such as /brain-expr/app_<id>/, and data/ sits next to index.html.
site_dir <- function(pathname) sub("app_[^/]*/.*$", "", sub("[^/]*$", "", pathname))

# Keeps the current selection when it is still offered, else the first choice.
keep_choice <- function(current, choices) if (isTRUE(current %in% choices)) current else choices[1]

parse_genes <- function(txt) {
  g <- unique(strsplit(trimws(txt), "[[:space:],;]+")[[1]])
  g[nzchar(g)]
}

fmt_int <- function(x) format(as.integer(x), big.mark = " ")

how_to_read <- function(...) div(class = "how-to-read", tags$b("Comment lire ce graphique. "), ...)

theme_app <- function() {
  theme_minimal(14) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 13, hjust = 0),
          legend.title = element_text(size = 12), legend.text = element_text(size = 11),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(colour = "grey35"))
}

css <- "
.how-to-read { background: #f4f8fb; border-left: 4px solid #2c6e9b; padding: 10px 14px;
               margin: 4px 0 14px; font-size: 0.95rem; }
.step { font-weight: 600; color: #2c6e9b; margin-bottom: 4px; }
.atlas-card { background: #f8f9fa; border-radius: 6px; padding: 8px 10px; font-size: 0.88rem; }
.missing { color: #b03a2e; font-size: 0.9rem; }
.btn-example { font-size: 0.8rem; padding: 2px 8px; margin: 2px 4px 2px 0; }
.guide h5 { color: #2c6e9b; margin-top: 1.2rem; }
"

ui <- page_sidebar(
  title = div("Expression de gènes dans le cerveau en développement",
              span(class = "fs-6 fw-normal ms-2 opacity-75", "souris et humain, cellule par cellule")),
  theme = bs_theme(version = 5, primary = "#2c6e9b", "navbar-bg" = "#1f4e6e",
                   base_font = font_collection("system-ui", "-apple-system", "Segoe UI", "Roboto", "sans-serif")),
  # Plots size themselves to the number of stages and cell types; a fillable
  # page would squeeze them to the window height.
  fillable = FALSE,
  tags$head(tags$style(HTML(css))),
  sidebar = sidebar(width = 340,
    div(class = "step", "1. Choisir un atlas"),
    selectInput("atlas", NULL, NULL),
    uiOutput("atlas_card"),
    div(class = "step mt-3", "2. Entrer vos gènes"),
    textAreaInput("genes_txt", NULL, EXAMPLES$rho_core, rows = 4,
                  placeholder = "ex. Trio Rac1 Cdc42"),
    div(class = "small text-muted mb-1",
        "Séparés par des espaces ou des virgules. Majuscules ou minuscules : peu importe. ",
        sprintf("Jusqu'à %d gènes.", MAX_GENES)),
    selectInput("preset", "Ou charger une liste prête", PRESETS),
    actionButton("go", "Afficher", class = "btn-primary w-100 mt-2"),
    uiOutput("gene_msg"),
    div(class = "step mt-3", "3. Affiner (facultatif)"),
    uiOutput("region_ui"),
    textInput("pos_gene", "Seulement les cellules qui expriment…", "",
              placeholder = "ex. Pvalb, Sst, Nkx2-1, Lhx6"),
    div(class = "small text-muted mb-2",
        "Un gène : l'app ne garde que les cellules où il est détecté. Appliqué avec « Afficher »."),
    selectInput("level", "Grouper les cellules par", NULL),
    checkboxGroupInput("stages", "Stades à afficher", NULL),
    checkboxInput("pool", "Regrouper les stades cochés en un seul", FALSE),
    accordion(open = FALSE,
      accordion_panel("Plus d'options",
        sliderInput("minpct", "Cacher les groupes où aucun gène n'atteint ce % de cellules",
                    0, 50, 5, post = " %"),
        textInput("ctfilter", "Ne montrer que les groupes contenant…", "",
                  placeholder = "ex. glia, neuron"),
        radioButtons("colour", "Couleur des points",
                     c("Niveau d'expression" = "mean_lognorm",
                       "Niveau relatif (chaque gène sur sa propre échelle)" = "z"))))
  ),
  navset_card_tab(id = "tabs",
    nav_panel("Vue d'ensemble",
      how_to_read("Chaque ligne est un groupe de cellules, chaque colonne un gène. ",
                  tags$b("Plus le point est gros"), ", plus la proportion de cellules du groupe qui expriment le gène est grande. ",
                  tags$b("Plus il est clair (jaune)"), ", plus le gène y est exprimé fort en moyenne. ",
                  "Un gros point jaune = gène exprimé par beaucoup de cellules, et fort. ",
                  "Les groupes qui se ressemblent sont placés côte à côte. ",
                  sprintf("Les groupes de moins de %d cellules sont cachés.", MIN_CELLS)),
      uiOutput("dot_ui")),
    nav_panel("Carte de chaleur",
      how_to_read("Une ligne par gène, une colonne par groupe de cellules, un panneau par stade. ",
                  tags$b("Plus la case est foncée"), ", plus la part des cellules du groupe qui expriment le gène est grande. ",
                  "Cette vue montre toute la liste (jusqu'à ", MAX_GENES, " gènes) ; ",
                  "pour les régulateurs des Rho GTPases, les gènes sont regroupés par famille."),
      uiOutput("heat_ui")),
    nav_panel("Au fil du développement",
      layout_columns(col_widths = c(4, 8),
        selectInput("time_gene", "Gène", NULL),
        radioButtons("time_metric", "Mesure", inline = TRUE,
                     c("% de cellules" = "pct_expressing", "Niveau moyen" = "mean_lognorm"))),
      how_to_read("Une ligne par groupe de cellules. Les 6 groupes où le gène est le plus exprimé ",
                  "sont en couleur et nommés ; les autres sont en gris. ",
                  "Une ligne qui monte = le gène est exprimé par de plus en plus de cellules au fil des stades. ",
                  "Un point manquant = moins de ", MIN_CELLS, " cellules de ce groupe à ce stade."),
      uiOutput("time_ui")),
    nav_panel("Coexpression",
      layout_columns(col_widths = c(6, 6),
        selectInput("cx_stage", "Stade", NULL),
        selectInput("cx_ct", "Groupe de cellules", NULL)),
      how_to_read("Pour chaque paire de gènes, dans le groupe choisi : les deux gènes sont-ils détectés ",
                  tags$b("dans les mêmes cellules"), " plus souvent que par hasard ? ",
                  tags$b("Rouge"), " = plus souvent (×2 = deux fois plus que le hasard), ",
                  tags$b("bleu"), " = moins souvent, blanc = comme le hasard. ",
                  "Les étoiles indiquent un écart statistiquement solide. ",
                  "Attention : le scRNA-seq rate une partie des ARN peu abondants, ",
                  "donc une case blanche ou bleue ne prouve pas que les gènes ne sont jamais ensemble."),
      plotOutput("coex", height = "560px", fill = FALSE),
      h6(class = "mt-3", "Détail des paires"),
      DTOutput("coex_tab")),
    nav_panel("Tableau",
      p(class = "text-muted", "Les chiffres derrière la vue d'ensemble, pour Excel ou R."),
      downloadButton("dl", "Télécharger (CSV)", class = "btn-outline-primary mb-2"),
      DTOutput("tab")),
    nav_panel("Guide", div(class = "guide", uiOutput("guide")))
  )
)

server <- function(input, output, session) {
  stores <- reactiveVal(NULL)

  # Opens every available store: a local directory in native R, otherwise the
  # data/ directory served next to the shinylive app.
  observe({
    root <- local_root()
    if (is.null(root)) {
      cd <- session$clientData
      req(cd$url_hostname)
      root <- paste0(cd$url_protocol, "//", cd$url_hostname,
                     if (nzchar(cd$url_port)) paste0(":", cd$url_port),
                     site_dir(cd$url_pathname), "data/")
    }
    src <- if (is_url(root)) paste0(root, ATLASES, "/") else file.path(root, ATLASES)
    s <- lapply(src, function(x) tryCatch(relabel_store(open_store(x)), error = function(e) NULL))
    names(s) <- ATLASES
    s <- Filter(Negate(is.null), s)
    validate(need(length(s) > 0, "Aucune donnée trouvée."))
    stores(s)
    updateSelectInput(session, "atlas",
                      choices = setNames(names(s), vapply(s, function(x) x$text$short, "")))
  }) |> bindEvent(session$clientData$url_hostname, once = TRUE)

  observeEvent(input$preset, {
    req(nzchar(input$preset))
    g <- preset_genes(input$preset, st()$meta$species)
    updateTextAreaInput(session, "genes_txt", value = paste(g, collapse = " "))
  })

  st <- reactive({ req(stores(), input$atlas); stores()[[input$atlas]] })

  observeEvent(st(), {
    s <- st()
    updateSelectInput(session, "level", choices = s$level_names, selected = s$level_names[1])
    stg <- unique(s$cells$stage[order(s$cells$stage_rank)])
    updateCheckboxGroupInput(session, "stages", choices = stg, selected = stg)
  })

  output$atlas_card <- renderUI({
    s <- st()
    div(class = "atlas-card", s$text$what,
        div(class = "text-muted mt-1",
            sprintf("%s cellules · %s gènes", fmt_int(s$meta$n_cells), fmt_int(s$meta$n_genes))))
  })

  # Store, grouping level and stages valid for the current atlas. After an atlas
  # switch the inputs still hold the previous atlas's values until the update
  # round-trips. With pooling, the checked stages become one stage.
  # Region choices for atlases that declare a region column.
  output$region_ui <- renderUI({
    s <- st()
    lv <- s$meta$filter_level
    if (is.null(lv) || is.na(lv)) return(NULL)
    r <- sort(unique(s$cells[[lv]]))
    selectizeInput("regions", "Régions", r, selected = NULL, multiple = TRUE,
                   options = list(placeholder = "Toutes les régions"))
  })

  # Cells expressing the "positive" gene, read once per click on Afficher.
  pos <- reactive({
    g <- trimws(input$pos_gene)
    if (!nzchar(g)) return(NULL)
    b <- read_genes(st(), g)
    validate(need(length(b) == 1, sprintf("« %s » est introuvable dans cet atlas.", g)))
    list(gene = names(b), cells = b[[1]]$cell)
  }) |> bindEvent(input$go, st(), ignoreNULL = FALSE)

  view <- reactive({
    s <- st()
    req(input$level %in% s$meta$levels)
    stg <- intersect(input$stages, s$cells$stage)
    req(length(stg) > 0)
    # Region and positive-gene filters drop cells by moving them out of every stage.
    drop <- rep(FALSE, nrow(s$cells))
    lv <- s$meta$filter_level
    if (!is.null(lv) && !is.na(lv) && length(input$regions)) drop <- drop | !s$cells[[lv]] %in% input$regions
    p <- pos()
    if (!is.null(p)) drop <- drop | !seq_len(nrow(s$cells)) %in% p$cells
    s$cells$stage[drop] <- NA
    if (isTRUE(input$pool) && length(stg) > 1) {
      hit <- s$cells$stage %in% stg
      s$cells$stage[hit] <- POOLED
      s$cells$stage_rank[hit] <- min(s$cells$stage_rank[hit])
      stg <- POOLED
    }
    list(store = s, level = input$level, stages = stg)
  })

  genes_req <- reactive(head(parse_genes(input$genes_txt), MAX_GENES)) |>
    bindEvent(input$go, ignoreNULL = FALSE)

  hits <- reactive(match_genes(st(), genes_req()))

  blocks_all <- reactive({
    g <- hits()$found
    validate(need(length(g) > 0, "Aucun des gènes demandés n'est dans cet atlas."))
    withProgress(message = sprintf("Lecture de %d gènes…", length(g)), value = 0.5, read_genes(st(), g))
  })
  # The dotplot, time course and coexpression show the first MAX_DOT genes.
  blocks <- reactive(head(blocks_all(), MAX_DOT))

  output$gene_msg <- renderUI({
    h <- hits()
    msg <- character()
    n <- length(parse_genes(input$genes_txt))
    if (n > MAX_GENES)
      msg <- c(msg, sprintf("Seuls les %d premiers gènes sont lus.", MAX_GENES))
    if (length(h$found) > MAX_DOT)
      msg <- c(msg, sprintf("Carte de chaleur : %d gènes. Vue d'ensemble, évolution et coexpression : les %d premiers.",
                            length(h$found), MAX_DOT))
    if (length(h$missing))
      msg <- c(msg, paste("Introuvable(s) dans cet atlas :", paste(h$missing, collapse = ", ")))
    if (length(msg)) div(class = "missing mt-2", HTML(paste(msg, collapse = "<br>")))
  })

  summ <- reactive({
    v <- view()
    d <- summarise_genes(v$store, blocks(), v$level, v$stages)
    d <- d[d$n_cells >= MIN_CELLS, ]
    if (nzchar(input$ctfilter))
      d <- d[grepl(input$ctfilter, d$celltype, ignore.case = TRUE, fixed = FALSE), ]
    keep <- names(which(tapply(d$pct_expressing, d$celltype, max) >= input$minpct))
    d <- d[d$celltype %in% keep, ]
    validate(need(nrow(d) > 0, "Aucun groupe ne passe les filtres : baissez le seuil dans « Plus d'options »."))
    d$z <- ave(d$mean_lognorm, d$gene, FUN = function(x) if (sd(x) > 0) (x - mean(x)) / sd(x) else 0)
    w <- tapply(d$mean_lognorm, list(d$celltype, paste(d$stage, d$gene)), mean)
    w[is.na(w)] <- 0
    ord <- if (nrow(w) > 2) rownames(w)[hclust(dist(w))$order] else rownames(w)
    d$celltype <- factor(d$celltype, levels = ord)
    d$gene <- factor(d$gene, levels = names(blocks()))
    d$stage <- factor(d$stage, levels = unique(d$stage[order(d$stage_rank)]))
    d
  })

  n_facets <- reactive(nlevels(summ()$stage))
  ncol_facets <- reactive(min(3L, n_facets()))

  output$dot <- renderPlot({
    d <- summ()
    rel <- input$colour == "z"
    ggplot(d, aes(gene, celltype)) +
      geom_point(aes(size = pct_expressing, colour = .data[[input$colour]])) +
      facet_wrap(~stage, ncol = ncol_facets()) +
      scale_size_area(max_size = 7, limits = c(0, 100), breaks = c(10, 25, 50, 75, 100),
                      labels = function(x) paste0(x, " %"),
                      name = "Cellules du groupe\nqui expriment le gène") +
      scale_colour_viridis_c(name = if (rel) "Niveau relatif\n(par gène)" else "Niveau moyen\nd'expression",
                             breaks = if (rel) c(-1.5, 0, 1.5) else waiver(),
                             labels = if (rel) c("plus bas", "moyen", "plus haut") else waiver()) +
      guides(size = guide_legend(order = 1), colour = guide_colourbar(order = 2)) +
      labs(x = NULL, y = NULL) + theme_app() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 12),
            axis.text.y = element_text(size = 11), panel.grid.major = element_line(colour = "grey92"))
  })
  # La hauteur dynamique passe par renderUI : renderPlot(height = function() ...)
  # boucle sans fin sous webR.
  output$dot_ui <- renderUI({
    rows <- ceiling(n_facets() / ncol_facets())
    plotOutput("dot", height = paste0(max(380, 110 + rows * (70 + 18 * nlevels(summ()$celltype))), "px"),
               fill = FALSE)
  })

  # --- carte de chaleur : toute la liste, groupée par famille Rho quand elle s'applique ---
  heat <- reactive({
    v <- view()
    d <- summarise_genes(v$store, blocks_all(), v$level, v$stages)
    d <- d[d$n_cells >= MIN_CELLS, ]
    if (nzchar(input$ctfilter))
      d <- d[grepl(input$ctfilter, d$celltype, ignore.case = TRUE), ]
    validate(need(nrow(d) > 0, "Aucun groupe d'au moins 20 cellules avec ces filtres."))
    fam <- RHO$family[match(toupper(d$gene), toupper(ifelse(nzchar(RHO$mouse), RHO$mouse, RHO$human)))]
    fam2 <- RHO$family[match(toupper(d$gene), toupper(RHO$human))]
    fam[is.na(fam)] <- fam2[is.na(fam)]
    d$family <- factor(ifelse(is.na(fam), "Votre liste", fam), c(unique(RHO$family), "Votre liste"))
    d$gene <- factor(d$gene, rev(names(blocks_all())))
    d$stage <- factor(d$stage, unique(d$stage[order(d$stage_rank)]))
    d$celltype <- factor(d$celltype)
    droplevels(d)
  })

  output$heat <- renderPlot({
    d <- heat()
    ggplot(d, aes(celltype, gene, fill = pct_expressing)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      facet_grid(family ~ stage, scales = "free", space = "free") +
      scale_fill_gradientn(colours = c("#f7fbff", "#9ecae1", "#3182bd", "#08306b"), limits = c(0, 100),
                           labels = function(x) paste0(x, " %"),
                           name = "Cellules du groupe\nqui expriment le gène") +
      labs(x = NULL, y = NULL) + theme_app() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
            axis.text.y = element_text(face = "italic", size = 9),
            strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
            panel.grid = element_blank(), panel.spacing = unit(4, "pt"),
            legend.position = "top", legend.justification = "left")
  })
  output$heat_ui <- renderUI({
    d <- heat()
    h <- max(420, 220 + 13 * nlevels(d$gene) + 14 * nlevels(d$family))
    w <- max(900, 260 + 16 * nlevels(d$celltype) * nlevels(d$stage))
    div(style = "overflow-x: auto;",
        plotOutput("heat", height = paste0(h, "px"), width = paste0(w, "px"), fill = FALSE))
  })

  # --- au fil du développement ---
  observeEvent(blocks(), {
    g <- names(blocks())
    updateSelectInput(session, "time_gene", choices = g, selected = keep_choice(isolate(input$time_gene), g))
  })

  output$time <- renderPlot({
    d <- summ()
    validate(need(nlevels(d$stage) > 1,
                  "Cet atlas (ou votre sélection) n'a qu'un stade : voir la vue d'ensemble."))
    req(input$time_gene %in% d$gene)
    d <- d[d$gene == input$time_gene, ]
    m <- input$time_metric
    top <- names(sort(tapply(d[[m]], as.character(d$celltype), max), decreasing = TRUE))[1:min(6, nlevels(droplevels(d$celltype)))]
    d$hl <- ifelse(as.character(d$celltype) %in% top, as.character(d$celltype), "autres")
    d$hl <- factor(d$hl, c(top, "autres"))
    last <- do.call(rbind, lapply(split(d[d$hl != "autres", ], d$celltype[d$hl != "autres"], drop = TRUE),
                                  function(x) x[which.max(x$stage_rank), ]))
    pal <- c(setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(top)], top),
             autres = "grey80")
    ggplot(d, aes(stage, .data[[m]], group = celltype, colour = hl)) +
      geom_line(data = d[d$hl == "autres", ], linewidth = 0.5) +
      geom_line(data = d[d$hl != "autres", ], linewidth = 1.2) +
      geom_point(data = d[d$hl != "autres", ], size = 2.2) +
      geom_text(data = last, aes(label = celltype), hjust = -0.08, size = 4, show.legend = FALSE) +
      scale_colour_manual(values = pal, guide = "none") +
      scale_x_discrete(expand = expansion(add = c(0.4, 3.5))) +
      scale_y_continuous(labels = if (m == "pct_expressing") function(x) paste0(x, " %") else waiver()) +
      labs(x = NULL, title = input$time_gene,
           y = if (m == "pct_expressing") "Cellules du groupe qui expriment le gène" else "Niveau moyen d'expression") +
      theme_app() + theme(plot.title = element_text(face = "bold.italic"),
                          axis.text.x = element_text(angle = 30, hjust = 1))
  })
  output$time_ui <- renderUI(plotOutput("time", height = "560px", fill = FALSE))

  # --- coexpression ---
  cx <- reactive({
    validate(need(length(blocks()) > 1, "Il faut au moins deux gènes."))
    v <- view()
    coexpression(v$store, blocks(), v$level, v$stages)
  })
  observeEvent(cx(), {
    s <- unique(cx()$stage[order(view()$store$cells$stage_rank[match(cx()$stage, view()$store$cells$stage)])])
    updateSelectInput(session, "cx_stage", choices = s, selected = keep_choice(isolate(input$cx_stage), s))
  })
  observeEvent(list(cx(), input$cx_stage), {
    d <- cx()[cx()$stage %in% input$cx_stage & cx()$n_cells >= MIN_CELLS, ]
    req(nrow(d) > 0)
    n <- tapply(d$n_cells, d$celltype, `[`, 1)
    ct <- names(sort(tapply(d$pct_both, d$celltype, max), decreasing = TRUE))
    updateSelectInput(session, "cx_ct", choices = setNames(ct, sprintf("%s (%s cellules)", ct, fmt_int(n[ct]))),
                      selected = keep_choice(isolate(input$cx_ct), ct))
  })
  cx_sel <- reactive({
    req(input$cx_stage, input$cx_ct)
    d <- cx()[cx()$stage == input$cx_stage & cx()$celltype == input$cx_ct, ]
    validate(need(nrow(d) > 0, "Choisissez un stade et un groupe."))
    d
  })

  stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

  output$coex <- renderPlot({
    d <- cx_sel()
    lv <- names(blocks())
    swap <- match(d$gene1, lv) > match(d$gene2, lv)
    a <- ifelse(swap, d$gene2, d$gene1); b <- ifelse(swap, d$gene1, d$gene2)
    ratio <- 2^d$log2_obs_exp
    lab <- ifelse(d$n_gene1 == 0 | d$n_gene2 == 0, "–",
                  paste0(ifelse(ratio >= 1, paste0("×", formatC(ratio, format = "f", digits = 1)),
                                paste0("÷", formatC(1 / ratio, format = "f", digits = 1))),
                         "\n", stars(d$p_adj)))
    p <- data.frame(a = factor(a, lv), b = factor(b, rev(lv)),
                    v = pmax(-2, pmin(2, d$log2_obs_exp)), lab = lab)
    ggplot(p, aes(a, b, fill = v)) + geom_tile(colour = "white", linewidth = 1) +
      geom_text(aes(label = lab), size = 3.6, lineheight = 0.85) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-2, 2),
                           breaks = c(-2, -1, 0, 1, 2), labels = c("÷4", "÷2", "hasard", "×2", "×4"),
                           name = "Ensemble dans\nles mêmes cellules,\npar rapport au hasard") +
      labs(x = NULL, y = NULL,
           title = sprintf("%s — %s", input$cx_ct, input$cx_stage),
           subtitle = sprintf("%s cellules dans ce groupe", fmt_int(d$n_cells[1])),
           caption = paste0("Étoiles : écart au hasard statistiquement solide (* p < 0,05 ; ** p < 0,01 ; *** p < 0,001 ; corrigé pour les tests multiples).\n",
                            "Sans étoile, l'écart peut être dû au hasard.  – : un des deux gènes n'est détecté dans aucune cellule du groupe.")) +
      coord_equal() + theme_app() +
      theme(axis.text = element_text(face = "italic", size = 12), panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1), plot.caption = element_text(hjust = 0, colour = "grey35"))
  })
  output$coex_tab <- renderDT({
    d <- cx_sel()
    d <- d[order(d$p_adj), ]
    out <- data.frame(`Gène A` = d$gene1, `Gène B` = d$gene2,
                      `Cellules avec A` = d$n_gene1, `Cellules avec B` = d$n_gene2,
                      `Avec les deux (observé)` = d$n_both,
                      `Avec les deux si hasard (attendu)` = round(d$expected_both, 1),
                      `Observé / attendu` = round(2^d$log2_obs_exp, 2),
                      `p corrigée` = signif(d$p_adj, 2), check.names = FALSE)
    datatable(out, rownames = FALSE, options = list(pageLength = 10, dom = "tp"))
  })

  # --- tableau ---
  tab <- reactive({
    d <- summ()
    data.frame(Atlas = st()$text$short, Stade = d$stage, `Groupe de cellules` = d$celltype,
               Gène = d$gene, `Nombre de cellules` = d$n_cells,
               `% qui expriment` = round(d$pct_expressing, 1),
               `Niveau moyen` = round(d$mean_lognorm, 3), check.names = FALSE)
  })
  output$tab <- renderDT(datatable(tab(), rownames = FALSE, filter = "top",
                                   options = list(pageLength = 25, scrollX = TRUE)))
  output$dl <- downloadHandler(
    filename = function() sprintf("expression_%s_%s.csv", input$atlas, input$level),
    content = function(file) utils::write.csv(tab(), file, row.names = FALSE))

  output$guide <- renderUI({
    s <- stores()
    tagList(
      h5("À quoi sert cette app"),
      p("Elle répond à la question : ", tags$em("dans quelles cellules du cerveau, dans quelle région et à quel stade, mes gènes sont-ils exprimés ?"),
        " Elle s'appuie sur des atlas publiés où l'ARN de chaque cellule a été séquencé une à une (scRNA-seq)."),
      h5("L'atlas « Souris · E10.5 → adulte, par région »"),
      p("Il réunit 8 études pour suivre le prosencéphale de E9 à l'adulte, région par région : ",
        "MGE, LGE, CGE, pallium, cortex somatosensoriel, hippocampe, striatum, amygdale."),
      tags$ul(
        tags$li(tags$b("Âges"), " : les cellules sont regroupées autour de E10.5, E12.5, E13.5, E15.5, P0, P10, P21, P40 et P100. ",
                "Le signe ≈ rappelle que les âges réels sont les plus proches disponibles. Pour les voir, groupez par « Âge exact »."),
        tags$li(tags$b("Types cellulaires harmonisés"), " : un même classement pour toutes les études, à partir de gènes marqueurs ",
                "(ex. Nkx2-1 pour les progeniteurs de la MGE, Lhx6 pour les interneurones issus de la MGE, Pvalb, Sst, Aqp4, Olig1, C1qa). ",
                "C'est une classification automatique, pas une annotation publiée : « Grouper par : Annotation publiée » montre celle des auteurs."),
        tags$li(tags$b("Cellules triées"), " : Allaway 2021 et Mayer 2018 ont trié leurs cellules sur un gène rapporteur (Dlx6, Lhx6). ",
                "Dans ces échantillons, les pourcentages décrivent la lignée triée, pas tout le tissu. « Grouper par : Tri des cellules » les isole."),
        tags$li(tags$b("Ce qui manque"), " : l'amygdale avant P21, le cortex somatosensoriel après P30, l'hippocampe avant E15.5, le striatum embryonnaire en dehors de la LGE. ",
                "P40 et P100 exacts n'existent dans aucun jeu public vérifié.")),
      if (!is.null(s$souris_dev)) p(class = "small text-muted", "Âges réels par fenêtre : ", s$souris_dev$meta$stage_note),
      h5("Régulateurs des Rho GTPases"),
      p("La liste « Régulateurs des Rho GTPases » (", nrow(RHO), " gènes) vient des familles officielles HGNC : Rho GTPases, GEF (familles Dbl et DOCK), ",
        "GAP, GDI et myosines IX, avec l'orthologue souris (MGI). La carte de chaleur les range par famille."),
      h5("Filtres"),
      tags$ul(
        tags$li(tags$b("Régions"), " : limite l'analyse aux régions choisies."),
        tags$li(tags$b("Seulement les cellules qui expriment…"), " : garde les cellules où ce gène est détecté (ex. Pvalb). ",
                "Le scRNA-seq ne détecte qu'une partie des ARN : ce filtre ne capture qu'une partie des vraies cellules positives. ",
                "Le gène du filtre lui-même apparaît alors à 100 %.")),
      p(tags$a(href = VIDEO_URL, target = "_blank", "▶ Voir la vidéo de démonstration (2 min)")),
      h5("En trois étapes"),
      tags$ol(
        tags$li(tags$b("Choisir un atlas"), " (souris ou humain, et la période de développement)."),
        tags$li(tags$b("Entrer ses gènes"), " puis cliquer sur ", tags$b("Afficher"),
                ". Les gènes introuvables dans l'atlas sont signalés en rouge."),
        tags$li(tags$b("Lire les onglets"), " : vue d'ensemble, évolution au fil des stades, coexpression, tableau.")),
      h5("Les quatre atlas"),
      tags$table(class = "table table-sm",
        tags$thead(tags$tr(tags$th("Atlas"), tags$th("Ce qu'il contient"), tags$th("Stades"), tags$th("Référence"))),
        tags$tbody(lapply(s, function(x) tags$tr(
          tags$td(tags$b(x$text$short)), tags$td(x$text$what), tags$td(x$text$stages),
          tags$td(class = "small", x$meta$citation))))),
      h5("Ce que montrent les chiffres"),
      tags$ul(
        tags$li(tags$b("% de cellules qui expriment le gène"), " : parmi les cellules du groupe, la part où au moins une molécule d'ARN du gène a été détectée. C'est la mesure la plus robuste."),
        tags$li(tags$b("Niveau moyen d'expression"), " : quantité d'ARN du gène, ramenée à la taille de chaque cellule, sur une échelle logarithmique, moyennée sur tout le groupe. 0 = jamais détecté ; 1 à 2 = expression nette ; au-delà de 3 = très forte."),
        tags$li(tags$b("Coexpression"), " : nombre de cellules où les deux gènes sont détectés ensemble, comparé à ce qu'on attendrait si les deux gènes étaient indépendants. Un test statistique (hypergéométrique, corrigé pour le grand nombre de paires testées) dit si l'écart est solide.")),
      h5("Précautions"),
      tags$ul(
        tags$li("Comparez les profils ", tags$b("au sein d'un même atlas"), ". D'un atlas à l'autre, les techniques et les noms des types cellulaires diffèrent : les valeurs absolues ne sont pas comparables."),
        tags$li("Le scRNA-seq ", tags$b("rate une partie des ARN peu abondants"), ". Un gène faiblement exprimé peut sembler absent. Une absence de détection n'est donc pas une preuve d'absence."),
        tags$li("Les stades souris et humains ne se correspondent qu'à peu près : E12.5–E13.5 chez la souris ≈ premier trimestre chez l'humain ; P21 ≈ petite enfance."),
        tags$li(sprintf("Les groupes de moins de %d cellules à un stade donné sont cachés : sur si peu de cellules, un pourcentage ne veut rien dire.", MIN_CELLS)),
        tags$li("Pour les atlas humains, au plus 200 cellules ont été tirées au hasard par stade et par type cellulaire, pour que l'app reste légère. Les pourcentages restent justes ; la coexpression perd un peu de puissance."),
        tags$li("La recherche ignore majuscules et minuscules (Scn1a trouve SCN1A chez l'humain), mais un gène dont le nom diffère entre souris et humain sera signalé introuvable."),
        tags$li("Pour confirmer un résultat : hybridation in situ (Allen Developing Mouse Brain Atlas), immunomarquage ou données spatiales.")),
      h5("Petit lexique"),
      tags$dl(
        tags$dt("scRNA-seq"), tags$dd("Séquençage de l'ARN cellule par cellule."),
        tags$dt("E12.5, P21"), tags$dd("Jour embryonnaire 12,5 ; 21 jours après la naissance (souris)."),
        tags$dt("Semaines après la conception (pcw)"), tags$dd("Âge du fœtus humain compté depuis la conception."),
        tags$dt("Groupe de cellules"), tags$dd("Type cellulaire, région disséquée ou étude d'origine, selon le choix « Grouper les cellules par »."))
    )
  })
}

shinyApp(ui, server)
