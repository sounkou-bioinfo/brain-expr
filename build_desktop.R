#!/usr/bin/env Rscript
# Construit l'app desktop (Electron) à partir de app/ avec shinyelectron.
# Stratégie "shinylive" : R tourne en WebAssembly dans l'app, l'utilisateur
# n'a rien à installer (ni R, ni Node). Les tables de app/data/ sont incluses.
#
# Usage : Rscript build_desktop.R            # plateforme courante
#         Rscript build_desktop.R win|mac|linux
#
# Prérequis machine de build : R >= 4.5, Node >= 22, npm >= 11.5
#   shinyelectron::sitrep_shinyelectron()  pour vérifier
#   shinyelectron::install_nodejs()        pour installer Node sans droits admin
#
# Notes :
# - Un mac ne se construit que sur macOS, Windows plus fiablement sur Windows :
#   le workflow .github/workflows/desktop.yml fait les trois en CI.
# - shinyelectron est un prototype (« pas recommandé en production ») : tester
#   l'exécutable avant de l'envoyer.
# - Les builds ne sont pas signés : avertissement SmartScreen sous Windows,
#   clic droit > Ouvrir la première fois sous macOS.

args <- commandArgs(trailingOnly = TRUE)
platform <- if (length(args)) args[1] else NULL

if (!file.exists("app/data/expression_summary.csv.gz"))
  message("!! app/data/ vide : l'app desktop affichera les données de DÉMONSTRATION. ",
          "Lancez run_pipeline.R d'abord.")

shinyelectron::export(
  appdir = "app",
  destdir = "build",
  app_name = "Brain Gene Expression",
  runtime_strategy = "shinylive",
  platform = platform,
  overwrite = TRUE
)
message("Exécutables : build/electron-app/dist/")
