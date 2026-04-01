rm(list = ls())
source("config.R")
set.seed(CONFIG$seed)

source("01_utils.R")
source("02_analysis.R")
source("03_ergm.R")
source("04_plotting.R")

dir.create(CONFIG$fig_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load data
dat <- load_data(CONFIG)

# 2. Network stats & clustering
Deg <- get_network_stat(dat$G7, "degree")
AA  <- expand_matrix(dat$YL, dat$A7)
clu <- reg.SSP(AA, K = CONFIG$n_clusters,
               nstart = CONFIG$nstart, tau = CONFIG$tau)
YL_clu <- dat$YL %>%
  dplyr::left_join(
    data.frame(idr = as.numeric(colnames(AA)), Clus = clu$cluster),
    by = "idr")

# 3. Outcome
Y <- YL_clu %>%
  dplyr::left_join(Deg, by = "idr") %>%
  dplyr::mutate(
    y7      = as.integer(L7 > 1),
    y9      = as.integer(L9 > 1),
    deg.out = ifelse(is.na(deg.out), 0, deg.out),
    Job     = as.integer(G512 == 1)) %>%
  dplyr::left_join(dat$X_Group, by = "idr")
Y$Group[is.na(Y$Group)] <- 0
Y$Job[is.na(Y$Job)]     <- 0

# 4. Neighbour predictors
mnei_lonely <- find_mutual_ties(Y, dat$G7, lonely = TRUE)
mnei_deg    <- find_avg_neighbor_degree(Y, dat$G7, K = CONFIG$K)
deg_nl      <- find_mutual_ties(Y, dat$G7, lonely = FALSE)

# 5. GLM / GLMER
model_tag <- ifelse(CONFIG$use_lmer, "GLMER", "GLM")
cat(sprintf("Fitting %s models...\n", model_tag))

sep  <- lapply(CONFIG$thresholds, function(i)
  fit_separated(Y, i, mnei_deg, deg_nl, mnei_lonely,
                K = CONFIG$K, tq = CONFIG$trim_q,
                use_lmer = CONFIG$use_lmer))

allc <- lapply(CONFIG$thresholds, function(i)
  fit_allcontacts(Y, i, mnei_deg,
                  tq = CONFIG$trim_q,
                  use_lmer = CONFIG$use_lmer))

# 6. Results
B  <- extract_betas(sep, allc)
or <- function(df, nm)
  exp(df %>% dplyr::select(!!nm := est, low = asymp.LCL, high = asymp.UCL)) %>% round(2)

cat("\n--- Separated ORs ---\n")
print(cbind(or(B$ls_d,"nL_lo"), or(B$ls_C,"L_lo"),
            or(B$su_d,"nL_hi"), or(B$su_C,"L_hi")))
cat("\n--- All-Contacts ORs ---\n")
print(cbind(or(B$deg_l,"all_lo"), or(B$deg_m,"all_hi")))

# 7. Figure
fig_name <- sprintf("Estimate_%s.png", model_tag)
png(file.path(CONFIG$fig_dir, fig_name), h = 6, w = 6, units = "in", res = 100)
plot_marginal(B)
dev.off()

# 8. ERGM
if (isTRUE(CONFIG$run_ergm)) {
  cat("\nFitting ERGM...\n")
  ergm_results <- run_ergm_analysis(Y, dat$A7)
}

cat("\nDone.\n")