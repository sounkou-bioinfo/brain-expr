# Runs the app's server logic on every atlas and renders every output, without a browser.
# Usage: BRAIN_EXPR_DATA=<stores_dir> Rscript tests/test_app.R app
suppressPackageStartupMessages(library(shiny))
app_dir <- commandArgs(TRUE)[1]
G <- new.env()
old <- setwd(app_dir); sys.source("app.R", G); setwd(old)

genes <- paste(G$EXAMPLES$epilepsy, "NotAGene")
for (atlas in G$ATLASES) {
  testServer(G$server, {
    session$setInputs(atlas = atlas, genes_txt = genes, go = 1)
    s <- st()
    for (lv in s$meta$levels) for (pool in c(FALSE, TRUE)) {
      session$setInputs(level = lv, stages = unique(s$cells$stage), pool = pool, minpct = 5,
                        ctfilter = "", colour = "mean_lognorm", time_metric = "pct_expressing")
      stopifnot(nrow(summ()) > 0, "NotAGene" %in% hits()$missing, !"Scn1a" %in% hits()$missing)
      session$setInputs(time_gene = names(blocks())[1])
      invisible(output$dot); invisible(output$dot_ui); invisible(output$tab)
      if (nlevels(summ()$stage) > 1) invisible(output$time)
      session$setInputs(cx_stage = cx()$stage[1])
      session$setInputs(cx_ct = cx()$celltype[cx()$stage == input$cx_stage][1])
      invisible(output$coex); invisible(output$coex_tab)
      session$setInputs(colour = "z", time_metric = "mean_lognorm")
      invisible(output$dot)
      cat(sprintf("OK %-14s %-14s pool=%-5s stages=%d groups=%d genes=%d missing=%s\n", atlas, lv, pool,
                  nlevels(summ()$stage), nlevels(summ()$celltype), length(blocks()),
                  paste(hits()$missing, collapse = ",")))
    }
    invisible(output$guide)
  })
}
