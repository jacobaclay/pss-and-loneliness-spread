### `config.R`

CONFIG <- list(
  seed       = 1232,
  data_dir   = "data/",
  fig_dir    = "output/",
  K          = 0.5,
  n_clusters = 8,
  nstart     = 50,
  tau        = 0.3,
  trim_q     = 0.005,
  thresholds = 11:7,
  use_lmer   = TRUE,   # TRUE = GLMER with (1|Clus); FALSE = standard GLM
  run_ergm   = TRUE
)