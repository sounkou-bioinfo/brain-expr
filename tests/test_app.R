# Runs the app's server logic on every atlas and renders every output, without a browser.
# Usage: BRAIN_EXPR_DATA=<stores_dir> Rscript tests/test_app.R app
suppressPackageStartupMessages(library(shiny))
app_dir <- commandArgs(TRUE)[1]
G <- new.env()
old <- setwd(app_dir); sys.source("app.R", G); setwd(old)

genes <- paste(G$EXAMPLES$rho_core, G$EXAMPLES$lineage, "NotAGene")
for (atlas in G$ATLASES) {
  testServer(G$server, {
    session$setInputs(atlas = atlas, genes_txt = genes, pos_gene = "", go = 1)
    s <- st()
    for (lv in s$meta$levels) for (pool in c(FALSE, TRUE)) {
      session$setInputs(level = lv, stages = unique(s$cells$stage), pool = pool, minpct = 5,
                        ctfilter = "", colour = "mean_lognorm", time_metric = "pct_expressing")
      stopifnot(nrow(summ()) > 0, "NotAGene" %in% hits()$missing, !"Trio" %in% hits()$missing)
      session$setInputs(time_gene = names(blocks())[1])
      invisible(output$dot); invisible(output$dot_ui); invisible(output$tab)
      invisible(output$heat); invisible(output$heat_ui)
      if (nlevels(summ()$stage) > 1) invisible(output$time)
      session$setInputs(cx_stage = cx()$stage[1])
      session$setInputs(cx_ct = cx()$celltype[cx()$stage == input$cx_stage][1])
      invisible(output$coex); invisible(output$coex_tab)
      cat(sprintf("OK %-14s %-10s pool=%-5s stages=%d groups=%d genes=%d\n", atlas, lv, pool,
                  nlevels(summ()$stage), nlevels(summ()$celltype), length(blocks())))
    }

    # Preset lists: the full Rho regulator list loads, the heatmap groups it by family.
    rho <- G$preset_genes("rho:all", s$meta$species)
    stopifnot(length(rho) == nrow(G$RHO))
    session$setInputs(genes_txt = paste(rho, collapse = " "), go = 2, level = s$meta$levels[1], pool = FALSE)
    stopifnot(length(blocks_all()) > 140, length(blocks()) == G$MAX_DOT,
              nlevels(heat()$family) >= 5)
    invisible(output$heat)
    cat(sprintf("OK %-14s Rho list: %d genes read, %d families in the heatmap\n", atlas,
                length(blocks_all()), nlevels(heat()$family)))

    # Positive-gene filter: every kept cell expresses the gene, so it shows 100 %.
    pg <- if (s$meta$species == "human") "GAD1" else "Gad1"
    session$setInputs(genes_txt = paste(pg, "Trio"), pos_gene = pg, go = 3)
    d <- summ()
    stopifnot(all(abs(d$pct_expressing[toupper(d$gene) == toupper(pg)] - 100) < 1e-9))
    cat(sprintf("OK %-14s filter %s+: %d groups\n", atlas, pg, nlevels(d$celltype)))

    # Region filter, for atlases that declare one.
    lv <- s$meta$filter_level
    if (!is.null(lv) && !is.na(lv)) {
      r <- c("MGE", "Cortex somatosensoriel")
      session$setInputs(pos_gene = "", regions = r, go = 4, level = lv)
      stopifnot(all(as.character(summ()$celltype) %in% r))
      cat(sprintf("OK %-14s region filter: %s\n", atlas, paste(levels(summ()$celltype), collapse = ", ")))
    }
    invisible(output$guide)
  })
}
