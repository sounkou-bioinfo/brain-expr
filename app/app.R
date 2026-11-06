# Explorateur d'expression : n'importe quelle liste de gènes x types cellulaires x stades,
# dans quatre atlas scRNA-seq du cerveau en développement (souris et humain).
#
# Les données sont des « stores » par atlas (format : store.R), lus gène par gène.
# - R natif : dossier local, BRAIN_EXPR_DATA ou ../stores à côté de app/.
# - shinylive / app desktop : dossier data/ servi à côté de l'app, lu par requêtes HTTP Range.

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(ggplot2); library(DT)
})
source("store.R", local = TRUE)

ATLASES <- c("lamanno2021", "zeisel2018", "braun2023", "velmeshev2023")
EXAMPLE_GENES <- "Pax6 Sox2 Hes5 Mki67 Eomes Neurod2 Tbr1 Gad2 Olig2 Aqp4"
MAX_GENES <- 30L

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

ui <- page_sidebar(
  title = "Expression de gènes dans le cerveau en développement",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  # Plots size themselves to the number of stages and cell types; a fillable
  # page would squeeze them to the window height.
  fillable = FALSE,
  sidebar = sidebar(width = 320,
    selectInput("atlas", "Atlas", NULL),
    textAreaInput("genes_txt", "Gènes (espaces, virgules ou retours à la ligne)",
                  EXAMPLE_GENES, rows = 4),
    actionButton("go", "Afficher", class = "btn-primary"),
    uiOutput("gene_msg"),
    selectInput("level", "Niveau d'annotation", NULL),
    checkboxGroupInput("stages", "Stades", NULL),
    sliderInput("minpct", "Masquer les types où aucun gène n'atteint (%)", 0, 50, 5),
    textInput("ctfilter", "Filtrer les types cellulaires (regex)", ""),
    radioButtons("colour", "Couleur", c("Moyenne log(CP10k+1)" = "mean_lognorm",
                                        "z-score par gène" = "z"))
  ),
  uiOutput("atlas_note"),
  navset_card_tab(
    nav_panel("Dotplot", uiOutput("dot_ui"),
              helpText("Taille = % de cellules où le gène est détecté ; couleur = niveau moyen. ",
                       "Comparer les profils relatifs, pas les valeurs absolues entre atlas.")),
    nav_panel("Dans le temps", uiOutput("time_ui"),
              helpText("% de cellules qui expriment chaque gène, par stade et par type cellulaire.")),
    nav_panel("Coexpression",
              layout_columns(col_widths = c(4, 4, 4),
                selectInput("cx_stage", "Stade", NULL),
                selectInput("cx_ct", "Type cellulaire", NULL),
                selectInput("cx_metric", "Mesure", c("log2 observé/attendu" = "log2_obs_exp",
                                                     "% cellules co-exprimant" = "pct_both",
                                                     "Jaccard" = "jaccard"))),
              div(style = "margin-top: 12px", plotOutput("coex", height = "520px", fill = FALSE)),
              helpText("Observé = cellules où les deux gènes sont détectés ; attendu = indépendance (pA × pB). ",
                       "Le dropout du scRNA-seq sous-estime la coexpression : une absence n'est pas une preuve."),
              DTOutput("coex_tab")),
    nav_panel("Table", downloadButton("dl", "Télécharger (CSV)"), DTOutput("tab")),
    nav_panel("À propos", uiOutput("about"))
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
    s <- lapply(src, function(x) tryCatch(open_store(x), error = function(e) NULL))
    names(s) <- ATLASES
    s <- Filter(Negate(is.null), s)
    validate(need(length(s) > 0, "Aucune donnée trouvée."))
    stores(s)
    updateSelectInput(session, "atlas",
                      choices = setNames(names(s), vapply(s, function(x) x$meta$label, "")))
  }) |> bindEvent(session$clientData$url_hostname, once = TRUE)

  st <- reactive({ req(stores(), input$atlas); stores()[[input$atlas]] })

  observeEvent(st(), {
    s <- st()
    updateSelectInput(session, "level", choices = s$meta$levels, selected = s$meta$levels[1])
    stg <- unique(s$cells[order(s$cells$stage_rank), "stage"])
    updateCheckboxGroupInput(session, "stages", choices = stg, selected = stg)
  })

  # Level and stages valid for the current atlas; after an atlas switch the
  # inputs still hold the previous atlas's values until the update round-trips.
  sel <- reactive({
    s <- st()
    req(input$level %in% s$meta$levels)
    stg <- intersect(input$stages, s$cells$stage)
    req(length(stg) > 0)
    list(level = input$level, stages = stg)
  })

  genes_req <- reactive(head(parse_genes(input$genes_txt), MAX_GENES)) |>
    bindEvent(input$go, ignoreNULL = FALSE)

  hits <- reactive(match_genes(st(), genes_req()))

  blocks <- reactive({
    g <- hits()$found
    validate(need(length(g) > 0, "Aucun des gènes demandés n'est dans cet atlas."))
    withProgress(message = "Lecture des gènes", value = 0.5, read_genes(st(), g))
  })

  output$gene_msg <- renderUI({
    h <- hits()
    msg <- character()
    if (length(parse_genes(input$genes_txt)) > MAX_GENES)
      msg <- c(msg, sprintf("Seuls les %d premiers gènes sont affichés.", MAX_GENES))
    if (length(h$missing))
      msg <- c(msg, paste("Absents de cet atlas :", paste(h$missing, collapse = ", ")))
    if (length(msg)) div(class = "small text-danger mt-2", HTML(paste(msg, collapse = "<br>")))
  })

  output$atlas_note <- renderUI({
    m <- st()$meta
    div(class = "small text-muted mb-2",
        sprintf("%s — %s cellules, %s gènes. Échantillonnage : %s.", m$citation,
                format(as.integer(m$n_cells), big.mark = " "), format(as.integer(m$n_genes), big.mark = " "),
                m$sampling))
  })

  summ <- reactive({
    d <- summarise_genes(st(), blocks(), sel()$level, sel()$stages)
    if (nzchar(input$ctfilter)) d <- d[grepl(input$ctfilter, d$celltype, ignore.case = TRUE), ]
    keep <- names(which(tapply(d$pct_expressing, d$celltype, max) >= input$minpct))
    d <- d[d$celltype %in% keep, ]
    validate(need(nrow(d) > 0, "Aucune ligne avec ces filtres."))
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
  ncol_facets <- reactive(min(4L, n_facets()))

  output$dot <- renderPlot({
    d <- summ()
    ggplot(d, aes(gene, celltype)) +
      geom_point(aes(size = pct_expressing, colour = .data[[input$colour]])) +
      facet_wrap(~stage, ncol = ncol_facets()) +
      scale_size_area(max_size = 6, name = "% cellules") +
      scale_colour_viridis_c(name = if (input$colour == "z") "z" else "moy. log(CP10k+1)") +
      labs(x = NULL, y = NULL) + theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            strip.text = element_text(face = "bold"))
  })
  # La hauteur dynamique passe par renderUI : renderPlot(height = function() ...)
  # boucle sans fin sous webR.
  output$dot_ui <- renderUI({
    rows <- ceiling(n_facets() / ncol_facets())
    plotOutput("dot", height = paste0(max(350, 80 + rows * (60 + 16 * nlevels(summ()$celltype))), "px"),
               fill = FALSE)
  })

  output$time <- renderPlot({
    d <- summ()
    validate(need(nlevels(d$stage) > 1, "Un seul stade dans cet atlas : voir le dotplot."))
    ggplot(d, aes(stage, pct_expressing, group = celltype, colour = celltype)) +
      geom_line() + geom_point(size = 1.5) +
      facet_wrap(~gene, ncol = 3, scales = "free_y") +
      labs(x = NULL, y = "% de cellules", colour = NULL) + theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
            strip.text = element_text(face = "bold.italic"))
  })
  output$time_ui <- renderUI(
    plotOutput("time", height = paste0(250 + 230 * ceiling(nlevels(summ()$gene) / 3), "px"), fill = FALSE))

  # --- coexpression ---
  cx <- reactive({
    validate(need(length(blocks()) > 1, "Il faut au moins deux gènes."))
    coexpression(st(), blocks(), sel()$level, sel()$stages)
  })
  observeEvent(cx(), {
    s <- unique(cx()$stage)
    updateSelectInput(session, "cx_stage", choices = s, selected = keep_choice(isolate(input$cx_stage), s))
  })
  observeEvent(list(cx(), input$cx_stage), {
    d <- cx()[cx()$stage == input$cx_stage, ]
    ct <- names(sort(tapply(d$pct_both, d$celltype, max), decreasing = TRUE))
    updateSelectInput(session, "cx_ct", choices = ct, selected = keep_choice(isolate(input$cx_ct), ct))
  })
  cx_sel <- reactive({
    req(input$cx_stage, input$cx_ct)
    cx()[cx()$stage == input$cx_stage & cx()$celltype == input$cx_ct, ]
  })
  output$coex <- renderPlot({
    d <- cx_sel()
    validate(need(nrow(d) > 0, "Pas de données de coexpression."))
    m <- input$cx_metric
    lv <- names(blocks())
    swap <- match(d$gene1, lv) > match(d$gene2, lv)
    a <- ifelse(swap, d$gene2, d$gene1); b <- ifelse(swap, d$gene1, d$gene2)
    p <- data.frame(a = factor(a, lv), b = factor(b, rev(lv)), v = d[[m]],
                    star = ifelse(d$p_adj < 0.001, "***", ifelse(d$p_adj < 0.01, "**",
                                  ifelse(d$p_adj < 0.05, "*", ""))))
    digits <- switch(m, pct_both = 0, jaccard = 2, 1)
    g <- ggplot(p, aes(a, b, fill = v)) + geom_tile(colour = "white") +
      geom_text(aes(label = paste0(round(v, digits), star)), size = 3.2) +
      labs(x = NULL, y = NULL, fill = m, caption = "* p_adj < 0.05  ** < 0.01  *** < 0.001",
           title = sprintf("%s — %s (n = %d cellules)", input$cx_stage, input$cx_ct, d$n_cells[1])) +
      coord_equal() + theme_minimal(12) +
      theme(axis.text = element_text(face = "italic"), axis.text.x = element_text(angle = 45, hjust = 1))
    if (m == "log2_obs_exp") g + scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c0392b")
    else g + scale_fill_viridis_c()
  })
  output$coex_tab <- renderDT({
    d <- cx_sel()
    d <- d[order(d$p_adj), c("gene1", "gene2", "n_cells", "n_gene1", "n_gene2", "n_both",
                             "pct_both", "expected_both", "log2_obs_exp", "jaccard", "p_adj")]
    datatable(d, rownames = FALSE, options = list(pageLength = 10)) |>
      formatRound(c("pct_both", "expected_both", "log2_obs_exp"), 1) |>
      formatRound("jaccard", 3) |> formatSignif("p_adj", 3)
  })

  # --- table ---
  tab <- reactive({
    d <- summ()
    data.frame(atlas = input$atlas, level = input$level, stage = d$stage, celltype = d$celltype,
               gene = d$gene, n_cells = d$n_cells, pct_expressing = round(d$pct_expressing, 1),
               mean_lognorm = round(d$mean_lognorm, 3))
  })
  output$tab <- renderDT(datatable(tab(), rownames = FALSE, filter = "top",
                                   options = list(pageLength = 25, scrollX = TRUE)))
  output$dl <- downloadHandler(
    filename = function() sprintf("expression_%s_%s.csv", input$atlas, input$level),
    content = function(file) utils::write.csv(tab(), file, row.names = FALSE))

  output$about <- renderUI({
    s <- stores()
    tagList(
      h5("Atlas"),
      tags$ul(lapply(s, function(x) tags$li(strong(x$meta$label), " — ", x$meta$citation,
                                             ". Échantillonnage : ", x$meta$sampling, "."))),
      h5("Méthode"),
      p("Chaque gène est lu à la demande dans la matrice de comptes. Normalisation : ",
        "log(1 + 10 000 × comptes / total de la cellule). Pour chaque gène × stade × type cellulaire : ",
        "nombre de cellules, % de cellules avec au moins un UMI, moyenne normalisée sur toutes les cellules du groupe. ",
        "Coexpression : test hypergéométrique unilatéral, correction de Benjamini-Hochberg sur toutes les paires et tous les groupes."),
      h5("Noms de gènes"),
      p("La recherche ignore la casse : Pax6 trouve PAX6 dans les atlas humains. ",
        "Ce n'est pas une table d'orthologues : un gène dont le nom diffère entre souris et humain sera signalé absent."),
      h5("Limites"),
      tags$ul(
        tags$li("Les atlas diffèrent par protocole, profondeur et annotation : comparer les profils relatifs au sein d'un atlas."),
        tags$li("Les stades souris et humains ne se correspondent qu'approximativement (E12.5-E13.5 ≈ 1er trimestre ; P21 ≈ enfance)."),
        tags$li("Les atlas humains sont sous-échantillonnés : les % et moyennes par groupe restent non biaisés, la puissance de la coexpression baisse."),
        tags$li("Validation recommandée : Allen Developing Mouse Brain Atlas (ISH) ou données spatiales.")))
  })
}

shinyApp(ui, server)
