# brain-expr — où et quand ~10 gènes sont exprimés dans le cerveau de souris

Pipeline R + app Shiny pour répondre à la question : pour une liste de gènes,
dans quels types cellulaires et à quels stades (E12.5, E13.5, P21) sont-ils
exprimés ou co-exprimés ?

| Stade | Source scRNA-seq | Fichier |
|---|---|---|
| E12.5, E13.5 | La Manno et al. 2021, *Nature* 596:92 (E7–E18, chaque jour) | `dev_all.loom` (~19 Go) |
| P21 (fenêtre P19–P23) | Zeisel et al. 2018, *Cell* 174:999 (P12–P30) | `l5_all.loom` |

Les deux fichiers se téléchargent depuis les pages *Downloads* de
mousebrain.org (`/development` et `/adolescent`). Placez-les dans `data/raw/`.

## Démarrage

```bash
# dépendances
Rscript -e 'install.packages(c("hdf5r","data.table","ggplot2","shiny","bslib","DT"))'

# 1. vérifier les noms d'attributs (âge, types cellulaires, tissu)
Rscript inspect_loom.R data/raw/dev_all.loom
Rscript inspect_loom.R data/raw/l5_all.loom

# 2. mettre la liste de gènes dans config.R (et ajuster les attributs si besoin)

# 3. extraire + résumer  (lit seulement les lignes des gènes demandés)
Rscript run_pipeline.R            # -> results/ et app/data/

# 4. explorer
Rscript -e 'shiny::runApp("app")'
```

Sans `app/data/`, l'app tourne sur des **données de démonstration aléatoires**
(bannière rouge), uniquement pour vérifier l'interface.

## Sorties (`results/`)

- `expression_summary.csv.gz` : dataset × niveau d'annotation × stade × type cellulaire × gène —
  `n_cells`, `pct_expressing` (% cellules ≥1 UMI), `mean_lognorm` (moy. log(1+CP10k)),
  `mean_lognorm_expressing`, `z_within_dataset`.
- `coexpression.csv.gz` : paires de gènes par groupe — `n_both`, `pct_both`, `expected_both`
  (pA·pB·n), `log2_obs_exp`, `jaccard`, `p_enrich` (hypergéométrique unilatéral), `p_adj` (BH).
- `cells_per_group.csv` : effectifs et âges réellement inclus.
- `dotplot_<dataset>_<niveau>.pdf` : figures statiques.

## App desktop (recommandé pour un utilisateur sans R)

`shinyelectron` emballe l'app dans Electron avec R compilé en WebAssembly
(stratégie *shinylive*) : un `.exe`, `.dmg` ou `.AppImage` qu'on ouvre en double-clic,
hors ligne, sans rien installer. Les tables de `app/data/` sont incluses.

```bash
Rscript run_pipeline.R            # remplir app/data/ avec les vraies tables
Rscript build_desktop.R           # -> build/electron-app/dist/
```

Pour Windows et macOS, le plus simple est la CI : pousser le dépôt sur GitHub avec
`app/data/` commité, puis lancer le workflow `desktop` (ou pousser un tag `v1.0`) ;
les exécutables apparaissent dans *Actions > Artifacts* (et dans la Release pour un tag).

Testé ici sur Linux : build ~190 Mo, l'app démarre en ~25 s, tous les onglets
fonctionnent. Deux pièges webR corrigés dans `app.R` (à ne pas réintroduire) :
`renderPlot(height = function() <reactive>)` et un `reactive` nommé `cur`
provoquent une récursion infinie sous webR alors qu'ils marchent en R natif ;
et `fread()` sur un `.gz` exige R.utils, d'où la lecture via `gzfile()`.
shinyelectron est encore un prototype : tester l'exécutable avant de l'envoyer.

## Déployer l'app en ligne

- **shinyapps.io** : voir `deploy_shinyapps.R` (compte gratuit, token, puis `Rscript deploy_shinyapps.R`).
- **Docker / serveur du labo** : `docker build -t brain-expr . && docker run -p 3838:3838 brain-expr`.

L'app ne contient que les tables résumées (quelques Mo), jamais les .loom.

## Limites à garder en tête

1. **Deux atlas différents** (protocoles, souches, profondeur) : comparer les profils
   *relatifs* au sein de chaque jeu, pas les niveaux absolus entre E13.5 et P21.
2. **Types cellulaires non équivalents** entre l'atlas embryonnaire et l'atlas adolescent.
3. **Dropout** : le scRNA-seq rate une partie des transcrits peu exprimés ; l'absence de
   coexpression n'est pas une preuve d'absence. À confirmer par ISH (Allen Developing Mouse
   Brain Atlas : E11.5, E13.5, E15.5, E18.5, P4, P14, P28) ou spatial (MOSTA, Stereo-seq,
   E9.5–E16.5, incluant E12.5 et E13.5).
4. P21 exact est peu échantillonné chez Zeisel : la fenêtre ±2 j est réglable
   (`stage_window` dans `config.R`) et les âges inclus sont tracés dans `cells_per_group.csv`.

## Tests

```bash
python3 tests/make_test_looms.py tests/fixtures      # looms synthétiques (pip install loompy)
Rscript run_pipeline.R tests/config_test.R /tmp/res  # pipeline de bout en bout
Rscript tests/test_app.R app                          # logique serveur Shiny
xvfb-run -a python3 tests/test_electron.py build/electron-app/dist/linux-unpacked/brain-gene-expression /tmp/shots
```
