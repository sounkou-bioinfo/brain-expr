# Explorateur d'expression : gènes x types cellulaires x stades (cerveau de souris)
# Lit les tables produites par run_pipeline.R dans app/data/.
# Sans données, l'app génère un jeu de DÉMONSTRATION aléatoire (bannière rouge).

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(data.table); library(ggplot2); library(DT)
})

data_dir <- "data"
f_sum  <- file.path(data_dir, "expression_summary.csv.gz")
f_coex <- file.path(data_dir, "coexpression.csv.gz")
f_prov <- file.path(data_dir, "provenance.txt")

make_demo <- function() {
  set.seed(42)
  genes <- c("Pax6", "Sox2", "Hes5", "Mki67", "Eomes", "Neurod2", "Tbr1", "Gad2", "Olig2", "Aqp4")
  spec <- list(
    lamanno2021 = list(stages = c("E12.5", "E13.5"),
                       ct = c("Radial glia", "Neuroblast", "Neuron", "Glioblast", "Vascular", "Immune")),
    zeisel2018  = list(stages = "P21",
                       ct = c("Neurons", "Astrocytes", "Oligodendrocytes", "OPC", "Microglia", "Vascular")))
  s <- rbindlist(lapply(names(spec), function(d) CJ(dataset = d, stage = spec[[d]]$stages,
                                                     celltype = spec[[d]]$ct, gene = genes)))
  s[, level := "Class"]
  s[, n_cells := sample(50:3000, 1), by = .(dataset, stage, celltype)]
  s[, pct_expressing := round(100 * rbeta(.N, 0.6, 2), 1)]
  s[, mean_lognorm := round(pct_expressing / 100 * runif(.N, 0.5, 2.5), 3)]
  s[, mean_lognorm_expressing := round(mean_lognorm / pmax(pct_expressing / 100, 0.01), 3)]
  s[, z_within_dataset := as.numeric(scale(mean_lognorm)), by = .(dataset, gene)]
  pr <- t(combn(genes, 2))
  cx <- unique(s[, .(dataset, level, stage, celltype, n_cells)])[, .(gene1 = pr[, 1], gene2 = pr[, 2]),
                                                                 by = .(dataset, level, stage, celltype, n_cells)]
  cx[s, p1 := i.pct_expressing / 100, on = .(dataset, stage, celltype, gene1 = gene)]
  cx[s, p2 := i.pct_expressing / 100, on = .(dataset, stage, celltype, gene2 = gene)]
  cx[, `:=`(n_gene1 = round(p1 * n_cells), n_gene2 = round(p2 * n_cells))]
  cx[, n_both := pmin(n_gene1, n_gene2, round(p1 * p2 * n_cells * runif(.N, 0.3, 3)))]
  cx[, `:=`(pct_both = 100 * n_both / n_cells, expected_both = n_gene1 * n_gene2 / n_cells)]
  cx[, log2_obs_exp := log2((n_both + 0.5) / (expected_both + 0.5))]
  cx[, jaccard := fifelse(n_gene1 + n_gene2 - n_both > 0, n_both / (n_gene1 + n_gene2 - n_both), 0)]
  cx[, p_enrich := phyper(n_both - 1, n_gene1, n_cells - n_gene1, n_gene2, lower.tail = FALSE)]
  cx[, p_adj := p.adjust(p_enrich, "BH"), by = dataset]
  cx[, c("p1", "p2") := NULL]
  list(sum = s, coex = cx)
}

# lecture .csv.gz sans R.utils (indisponible dans webR/shinylive)
read_gz <- function(f) { con <- gzfile(f, "rt"); on.exit(close(con)); fread(text = readLines(con)) }
DEMO <- !file.exists(f_sum)
if (DEMO) { d <- make_demo(); SUM <- d$sum; COEX <- d$coex } else {
  SUM <- read_gz(f_sum); COEX <- if (file.exists(f_coex)) read_gz(f_coex) else NULL
}
PROV <- if (file.exists(f_prov)) readLines(f_prov) else character()
TIME_H <- max(400, 140 + 16 * SUM[, uniqueN(celltype)])
ALL_GENES <- unique(SUM$gene)
DS_LABEL <- c(lamanno2021 = "La Manno 2021 — embryonnaire",
              zeisel2018  = "Zeisel 2018 — adolescent")
ds_choices <- setNames(unique(SUM$dataset),
                       ifelse(is.na(DS_LABEL[unique(SUM$dataset)]), unique(SUM$dataset),
                              DS_LABEL[unique(SUM$dataset)]))

ui <- page_sidebar(
  title = "Expression de gènes dans le cerveau de souris en développement",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(width = 300,
    selectInput("ds", "Jeu de données", ds_choices),
    selectInput("level", "Niveau d'annotation", NULL),
    checkboxGroupInput("stages", "Stades", NULL),
    selectizeInput("genes", "Gènes", ALL_GENES, selected = ALL_GENES, multiple = TRUE,
                   options = list(plugins = list("remove_button"))),
    sliderInput("minpct", "Masquer les types où aucun gène n'atteint (%)", 0, 50, 5),
    textInput("ctfilter", "Filtrer les types cellulaires (regex)", ""),
    radioButtons("colour", "Couleur", c("Moyenne log(CP10k+1)" = "mean_lognorm",
                                        "z-score par gène" = "z_within_dataset"), inline = FALSE)
  ),
  if (DEMO) div(class = "alert alert-danger", strong("DONNÉES DE DÉMONSTRATION ALÉATOIRES."),
                " Lancez run_pipeline.R pour générer app/data/ à partir des vrais fichiers .loom."),
  navset_card_tab(
    nav_panel("Dotplot", uiOutput("dot_ui"),
              helpText("Taille = % de cellules où le gène est détecté ; couleur = niveau moyen. ",
                       "Comparer les profils relatifs, pas les valeurs absolues entre jeux de données.")),
    nav_panel("Dans le temps", plotOutput("time", height = paste0(TIME_H, "px")),
              helpText("Évolution par stade, tous jeux de données confondus ; les types cellulaires ",
                       "sont ceux de chaque atlas et ne correspondent pas un à un entre embryon et P21.")),
    nav_panel("Coexpression",
              layout_columns(col_widths = c(4, 4, 4),
                selectInput("cx_stage", "Stade", NULL),
                selectInput("cx_ct", "Type cellulaire", NULL),
                selectInput("cx_metric", "Mesure", c("log2 observé/attendu" = "log2_obs_exp",
                                                     "% cellules co-exprimant" = "pct_both",
                                                     "Jaccard" = "jaccard"))),
              div(style = "margin-top: 12px", plotOutput("coex", height = "520px")),
              helpText("Observé = cellules où les deux gènes sont détectés ; attendu = indépendance (pA × pB). ",
                       "Le dropout du scRNA-seq sous-estime la coexpression : une absence n'est pas une preuve."),
              DTOutput("coex_tab")),
    nav_panel("Table", downloadButton("dl", "Télécharger (CSV)"), DTOutput("tab")),
    nav_panel("À propos",
      h5("Sources"),
      tags$ul(
        tags$li("La Manno et al. 2021, Nature 596:92 — atlas scRNA-seq du cerveau embryonnaire, E7–E18 (mousebrain.org/development)."),
        tags$li("Zeisel et al. 2018, Cell 174:999 — système nerveux adolescent, P12–P30 (mousebrain.org/adolescent)."),
        tags$li("Validation recommandée : Allen Developing Mouse Brain Atlas (ISH, E11.5–P28) et MOSTA (Stereo-seq, E9.5–E16.5).")),
      h5("Méthode"),
      p("Normalisation log(1 + 10 000 × comptes / total de la cellule). Pour chaque gène × type cellulaire × stade : ",
        "nombre de cellules, % de cellules avec ≥1 UMI, moyenne normalisée. Coexpression : test hypergéométrique unilatéral, BH par jeu de données."),
      h5("Provenance"), pre(if (length(PROV)) paste(PROV, collapse = "\n") else "démonstration"))
  )
)

server <- function(input, output, session) {
  observeEvent(input$ds, {
    lv <- unique(SUM[dataset == input$ds, level])
    updateSelectInput(session, "level", choices = lv, selected = lv[1])
    st <- sort(unique(SUM[dataset == input$ds, stage]))
    updateCheckboxGroupInput(session, "stages", choices = st, selected = st)
  })
  sel_data <- reactive({
    req(input$level, input$stages, input$genes)
    d <- SUM[dataset == input$ds & level == input$level & stage %in% input$stages & gene %in% input$genes]
    if (nzchar(input$ctfilter)) d <- d[grepl(input$ctfilter, celltype, ignore.case = TRUE)]
    keep <- d[, .(mx = max(pct_expressing)), by = celltype][mx >= input$minpct, celltype]
    d <- d[celltype %in% keep]
    validate(need(nrow(d) > 0, "Aucune ligne avec ces filtres."))
    # ordonner les types par regroupement hiérarchique du profil
    w <- dcast(d, celltype ~ stage + gene, value.var = "mean_lognorm", fill = 0)
    ord <- if (nrow(w) > 2) w$celltype[hclust(dist(as.matrix(w[, -1])))$order] else w$celltype
    d[, celltype := factor(celltype, levels = ord)]
    d[, gene := factor(gene, levels = intersect(input$genes, unique(gene)))]
    d
  })
  dot_h <- reactive(max(350, 120 + 18 * uniqueN(sel_data()$celltype)))
  output$dot <- renderPlot({
    d <- sel_data()
    ggplot(d, aes(gene, celltype)) +
      geom_point(aes(size = pct_expressing, colour = .data[[input$colour]])) +
      facet_wrap(~stage, nrow = 1) +
      scale_size_area(max_size = 7, name = "% cellules") +
      scale_colour_viridis_c(name = if (input$colour == "mean_lognorm") "moy. log(CP10k+1)" else "z") +
      labs(x = NULL, y = NULL) + theme_bw(13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            strip.text = element_text(face = "bold"))
  })
  # hauteur dynamique via renderUI : renderPlot(height = function() <reactive>)
  # part en récursion infinie sous webR (shinylive / shinyelectron)
  output$dot_ui <- renderUI(plotOutput("dot", height = paste0(dot_h(), "px")))

  output$time <- renderPlot({
    req(input$genes)
    # un niveau par jeu : celui choisi pour le jeu courant, le premier sinon
    lv <- SUM[, .(lv = unique(level)[1]), by = dataset]
    lv[dataset == input$ds, lv := input$level]
    d <- SUM[gene %in% input$genes][lv, on = .(dataset, level = lv), nomatch = 0]
    if (nzchar(input$ctfilter)) d <- d[grepl(input$ctfilter, celltype, ignore.case = TRUE)]
    d[, stage_day := as.numeric(sub("^[EP]", "", stage)) + ifelse(grepl("^P", stage), 19.5, 0)]
    d[, stage := factor(stage, levels = unique(stage[order(stage_day)]))]
    d[, gene := factor(gene, levels = input$genes)]
    ggplot(d, aes(stage, celltype)) +
      geom_tile(aes(fill = mean_lognorm), colour = "white") +
      geom_text(aes(label = round(pct_expressing),
                    colour = mean_lognorm > 0.6 * max(mean_lognorm)), size = 2.6, show.legend = FALSE) +
      scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
      facet_grid(dataset ~ gene, scales = "free_y", space = "free_y") +
      scale_fill_viridis_c(name = "moy. log(CP10k+1)", option = "magma", direction = -1) +
      labs(x = NULL, y = NULL, caption = "chiffre = % de cellules exprimant") + theme_bw(11) +
      theme(strip.text.x = element_text(face = "italic"), axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # --- coexpression ---
  cx_cur <- reactive({
    req(!is.null(COEX), input$level)
    COEX[dataset == input$ds & level == input$level &
         gene1 %in% input$genes & gene2 %in% input$genes]
  })
  observe({
    d <- cx_cur(); st <- sort(unique(d$stage))
    updateSelectInput(session, "cx_stage", choices = st, selected = or_default(isolate(input$cx_stage), st[1]))
  })
  observe({
    req(input$cx_stage)
    d <- cx_cur()[stage == input$cx_stage]
    ct <- d[, .(m = max(pct_both)), by = celltype][order(-m), celltype]
    updateSelectInput(session, "cx_ct", choices = ct, selected = ct[1])
  })
  cx_sel <- reactive({
    req(input$cx_stage, input$cx_ct)
    cx_cur()[stage == input$cx_stage & celltype == input$cx_ct]
  })
  output$coex <- renderPlot({
    d <- cx_sel(); validate(need(nrow(d) > 0, "Pas de données de coexpression."))
    m <- input$cx_metric
    full <- d[, .(a = gene1, b = gene2, v = get(m), padj = p_adj)]
    lv <- intersect(input$genes, unique(c(full$a, full$b)))
    # triangle : a avant b dans l'ordre choisi
    swap <- match(full$a, lv) > match(full$b, lv)
    full[swap, c("a", "b") := .(b, a)]
    full[, `:=`(a = factor(a, lv), b = factor(b, rev(lv)))]
    full[, star := fifelse(padj < 0.001, "***", fifelse(padj < 0.01, "**", fifelse(padj < 0.05, "*", "")))]
    lab_fmt <- if (m == "pct_both") 0 else if (m == "jaccard") 2 else 1
    p <- ggplot(full, aes(a, b, fill = v)) + geom_tile(colour = "white") +
      geom_text(aes(label = paste0(round(v, lab_fmt), star)), size = 3.2) +
      labs(x = NULL, y = NULL, fill = names(which(c("log2 O/E" = "log2_obs_exp", "% both" = "pct_both",
                                                      "Jaccard" = "jaccard") == m)),
           title = sprintf("%s — %s (n = %s cellules)", input$cx_stage, input$cx_ct, d$n_cells[1])) +
      coord_equal() + theme_minimal(12) + labs(caption = "* p_adj<0.05  ** <0.01  *** <0.001") +
      theme(axis.text = element_text(face = "italic"), axis.text.x = element_text(angle = 45, hjust = 1))
    if (m == "log2_obs_exp") p + scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c0392b")
    else p + scale_fill_viridis_c()
  })
  output$coex_tab <- renderDT({
    d <- cx_sel()[order(p_adj), .(gene1, gene2, n_cells, n_gene1, n_gene2, n_both,
                                  pct_both = round(pct_both, 1), expected = round(expected_both, 1),
                                  log2_OE = round(log2_obs_exp, 2), jaccard = round(jaccard, 3),
                                  p_adj = signif(p_adj, 3))]
    datatable(d, rownames = FALSE, options = list(pageLength = 10))
  })

  # --- table ---
  tab <- reactive(sel_data()[, .(dataset, level, stage, celltype, gene, n_cells,
                            pct_expressing = round(pct_expressing, 1),
                            mean_lognorm = round(mean_lognorm, 3),
                            mean_lognorm_expressing = round(mean_lognorm_expressing, 3))])
  output$tab <- renderDT(datatable(tab(), rownames = FALSE, filter = "top",
                                   options = list(pageLength = 25, scrollX = TRUE)))
  output$dl <- downloadHandler(
    filename = function() sprintf("expression_%s_%s.csv", input$ds, input$level),
    content = function(file) fwrite(tab(), file))
}

or_default <- function(a, b) if (is.null(a) || !length(a) || !nzchar(a[1])) b else a
shinyApp(ui, server)
