# Runs the app's server logic on every atlas and renders every output, without a browser.
# Usage: BRAIN_EXPR_DATA=<stores_dir> Rscript tests/test_app.R app
suppressPackageStartupMessages(library(shiny))
app_dir <- commandArgs(TRUE)[1]
G <- new.env()
old <- setwd(app_dir); sys.source("app.R", G); setwd(old)

genes <- "Pax6 Sox2 Hes5 Mki67 Eomes Neurod2 Tbr1 Gad2 Olig2 Aqp4 NotAGene"
for (atlas in G$ATLASES) {
  testServer(G$server, {
    session$setInputs(atlas = atlas, genes_txt = genes, go = 1)
    s <- st()
    for (lv in s$meta$levels) {
      stg <- unique(s$cells$stage)
      session$setInputs(level = lv, stages = stg, minpct = 5, ctfilter = "",
                        colour = "mean_lognorm", cx_metric = "log2_obs_exp")
      stopifnot(nrow(summ()) > 0, "NotAGene" %in% hits()$missing, !"Pax6" %in% hits()$missing)
      invisible(output$dot); invisible(output$dot_ui); invisible(output$time_ui); invisible(output$tab)
      session$setInputs(cx_stage = cx()$stage[1])
      session$setInputs(cx_ct = cx()$celltype[cx()$stage == input$cx_stage][1])
      invisible(output$coex); invisible(output$coex_tab)
      session$setInputs(colour = "z", cx_metric = "jaccard")
      invisible(output$dot); invisible(output$coex)
      cat(sprintf("OK %-14s %-14s genes=%d rows=%d types=%d pairs=%d missing=%s\n", atlas, lv,
                  length(blocks()), nrow(summ()), nlevels(summ()$celltype), nrow(cx_sel()),
                  paste(hits()$missing, collapse = ",")))
    }
  })
}
