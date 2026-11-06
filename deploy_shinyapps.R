# Déploiement sur shinyapps.io (gratuit jusqu'à 5 apps / 25 h actives par mois)
# 1. Créer un compte sur https://www.shinyapps.io
# 2. Account > Tokens > Show : copier name / token / secret ci-dessous (une seule fois)
# 3. Avoir lancé run_pipeline.R pour que app/data/ contienne les vraies tables
# 4. Rscript deploy_shinyapps.R

if (!requireNamespace("rsconnect", quietly = TRUE)) install.packages("rsconnect")
# rsconnect::setAccountInfo(name = "<compte>", token = "<token>", secret = "<secret>")

stopifnot("app/data/expression_summary.csv.gz absent : lancez d'abord run_pipeline.R" =
            file.exists("app/data/expression_summary.csv.gz"))
rsconnect::deployApp(appDir = "app", appName = "brain-gene-expression",
                     appFiles = c("app.R", file.path("data", list.files("app/data"))),
                     forceUpdate = TRUE)
