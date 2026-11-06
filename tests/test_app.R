# Exécute la logique serveur de l'app et rend chaque sortie (sans navigateur)
suppressPackageStartupMessages(library(shiny))
app_dir <- commandArgs(TRUE)[1]
G <- new.env(); old <- setwd(app_dir); sys.source("app.R", G); setwd(old)
SUM <- G$SUM
for (DS in unique(SUM$dataset)) for (LV in unique(SUM[dataset == DS, level])) {
  testServer(G$server, {
    session$setInputs(ds = DS)
    session$setInputs(level = LV, stages = unique(SUM[dataset == DS, stage]), genes = G$ALL_GENES,
                      minpct = 5, ctfilter = "", colour = "mean_lognorm", cx_metric = "log2_obs_exp")
    stopifnot(nrow(sel_data()) > 0)
    invisible(output$dot); invisible(output$time)
    session$setInputs(cx_stage = sort(unique(cx_cur()$stage))[1])
    session$setInputs(cx_ct = as.character(cx_cur()[stage == input$cx_stage, celltype][1]))
    invisible(output$coex); invisible(output$coex_tab); invisible(output$tab)
    session$setInputs(colour = "z_within_dataset", cx_metric = "pct_both")
    invisible(output$dot); invisible(output$coex)
    cat(sprintf("OK demo=%s %s/%s lignes=%d types=%d\n", G$DEMO, DS, LV, nrow(sel_data()), uniqueN(sel_data()$celltype)))
  })
}
