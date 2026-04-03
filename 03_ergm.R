source('config.R')

library(network)
library(ergm)
library(dplyr)

# ---------------------------------------------------------------------------
# Internal helper: subset Y and A to common IDs
# ---------------------------------------------------------------------------
.ergm_subset <- function(Y, A) {
  common <- intersect(as.character(Y$idr), colnames(A))
  Y_sub  <- Y %>%
    dplyr::filter(idr %in% as.numeric(common)) %>%
    dplyr::arrange(match(as.character(idr), common))
  A_sub  <- A[common, common]
  list(Y_sub = Y_sub, A_sub = A_sub, common = common)
}

# ---------------------------------------------------------------------------
# Build network object with vertex attributes
# ---------------------------------------------------------------------------
build_ergm_network <- function(Y, A) {
  s   <- .ergm_subset(Y, A)
  net <- network(s$A_sub, directed = TRUE, loops = FALSE,
                 matrix.type = "adjacency")
  net %v% "AGE"    <- s$Y_sub$AGE
  net %v% "SEX"    <- s$Y_sub$SEX
  net %v% "edu"    <- s$Y_sub$edu
  net %v% "spouse" <- s$Y_sub$spouse
  net %v% "D"      <- s$Y_sub$D
  net %v% "Job"      <- s$Y_sub$Job
  net %v% "Group"      <- s$Y_sub$Group
  net %v% "deg"    <- s$Y_sub$deg.out
  cat(sprintf("Network: %d nodes, %d edges\n",
              network.size(net), network.edgecount(net)))
  net
}

# ---------------------------------------------------------------------------
# Compute support-based edge-covariate matrices
# ---------------------------------------------------------------------------
build_support_matrices <- function(Y, A, N = 11) {
  s    <- .ergm_subset(Y, A)
  Supp <- s$Y_sub$B.love + s$Y_sub$B.listen + s$Y_sub$B.support
  S    <- as.integer(Supp >= N)
  S[is.na(S)] <- 0
  LS   <- 1L - S
  
  LS_mat <- outer(LS, LS, FUN = function(x, y) (x + y) > 0)
  S_mat  <- outer(S,  S,  FUN = function(x, y) (x + y) > 0)
  
  list(LS             = LS,
       S              = S,
       S_interaction  = s$A_sub * S_mat,
       LS_interaction = s$A_sub * LS_mat)
}

# ---------------------------------------------------------------------------
# Fit ERGM (direct model, no y7 as predictor)
# ---------------------------------------------------------------------------
fit_ergm_direct <- function(Y, A, N = 11, seed = 1232) {
  s    <- .ergm_subset(Y, A)
  net  <- build_ergm_network(Y, A)
  
  Supp <- s$Y_sub$B.love + s$Y_sub$B.listen + s$Y_sub$B.support
  S    <- as.integer(Supp >= N);  S[is.na(S)]   <- 0
  LS   <- 1L - S
  
  net %v% "LS" <- LS
  net %v% "S"  <- S
  
  keep  <- which(!is.na(LS))
  G_sub <- network::get.inducedSubgraph(net, v = keep)
  
  LS_sub <- network::get.vertex.attribute(G_sub, "LS")
  S_sub  <- network::get.vertex.attribute(G_sub, "S")
  Con    <- s$A_sub[keep, keep]
  
  LS_mat <- outer(LS_sub, LS_sub, FUN = function(x, y) (x + y) > 0)
  S_mat  <- outer(S_sub,  S_sub,  FUN = function(x, y) (x + y) > 0)
  
  ctrl <- control.ergm(
    seed            = seed,
    MCMC.burnin     = 50000,
    MCMC.interval   = 2000,
    MCMC.samplesize = 5000
  )
  
  ergm(
    G_sub ~ edges +
      nodecov("AGE")    +
      nodecov("SEX")    +
      nodecov("edu")    +
      nodecov("spouse") +
      nodecov("D")      +
      nodecov("Job") +
      nodecov("Group")      +
      edgecov(Con * S_mat)  +
      edgecov(Con * LS_mat) +
      nodecov("deg")    +
      nodematch("LS")   +
      nodecov("deg"):nodematch("LS"),
    control = ctrl
  )
}

# ---------------------------------------------------------------------------
# Tidy coefficient table
# ---------------------------------------------------------------------------
tidy_ergm <- function(fit) {
  cf <- summary(fit)$coefficients
  data.frame(
    term    = rownames(cf),
    est     = round(cf[, "Estimate"],    3),
    SE      = round(cf[, "Std. Error"],  3),
    p_val   = signif(cf[, "Pr(>|z|)"],  3),
    OR      = round(exp(cf[, "Estimate"]), 3),
    OR_low  = round(exp(cf[, "Estimate"] - 1.96 * cf[, "Std. Error"]), 3),
    OR_high = round(exp(cf[, "Estimate"] + 1.96 * cf[, "Std. Error"]), 3),
    row.names = NULL
  )
}

# ---------------------------------------------------------------------------
# GOF diagnostics
# ---------------------------------------------------------------------------
ergm_gof_check <- function(fit, fig_path = NULL, nsim = 100) {
  gof_res <- gof(fit, nsim = nsim)
  if (!is.null(fig_path)) {
    png(fig_path, height = 8, width = 10, units = "in", res = 150)
    par(mfrow = c(2, 2))
    plot(gof_res)
    dev.off()
  }
  invisible(gof_res)
}

# ---------------------------------------------------------------------------
# Sweep across support thresholds
# ---------------------------------------------------------------------------
ergm_sweep <- function(Y, A, thresholds = 11:7, seed = 1232) {
  lapply(setNames(thresholds, as.character(thresholds)), function(thr) {
    cat(sprintf("  ERGM threshold = %d\n", thr))
    fit <- tryCatch(
      fit_ergm_direct(Y, A, N = thr, seed = seed),
      error = function(e) {
        cat(sprintf("    failed: %s\n", e$message)); NULL
      }
    )
    if (is.null(fit)) return(NULL)
    list(threshold = thr, fit = fit, tidy = tidy_ergm(fit))
  })
}

# ---------------------------------------------------------------------------
# Run full ERGM analysis
# ---------------------------------------------------------------------------
run_ergm_analysis <- function(Y, A, cfg = CONFIG) {
  dir.create(cfg$fig_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n--- Fitting ERGM across thresholds ---\n")
  
  sweep_res <- ergm_sweep(Y, A,
                          thresholds = cfg$thresholds,
                          seed       = cfg$seed)
  sweep_res <- Filter(Negate(is.null), sweep_res)
  
  for (nm in names(sweep_res)) {
    cat(sprintf("\nThreshold = %s\n", nm))
    print(sweep_res[[nm]]$tidy)
  }
  
  mid <- as.character(cfg$thresholds[ceiling(length(cfg$thresholds) / 2)])
  if (mid %in% names(sweep_res)) {
    ergm_gof_check(sweep_res[[mid]]$fit,
                   file.path(cfg$fig_dir, "ergm_gof.png"))
  }
  
  invisible(sweep_res)
}
