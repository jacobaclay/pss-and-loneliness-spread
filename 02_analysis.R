library(lme4)
library(emmeans)
library(dplyr)

# Load and prepare data
# ===================================================================
load_data <- function(cfg) {
  e <- new.env()
  
  load(file.path(cfg$data_dir, "X_BS.RData"),    envir = e)
  X_Group <- e$X %>% dplyr::select(idr, Group)
  
  load(file.path(cfg$data_dir, "y_7_9.RData"),    envir = e)
  e$y_7_9$L7[is.nan(e$y_7_9$L7)] <- NA
  e$y_7_9$L9[is.nan(e$y_7_9$L9)] <- NA
  
  load(file.path(cfg$data_dir, "A7.RData"),        envir = e)
  
  mutu <- intersect(colnames(e$G0), e$X$idr)
  Y_mu <- e$X %>% dplyr::filter(idr %in% as.numeric(mutu))
  YL   <- e$X %>% dplyr::left_join(e$y_7_9, by = "idr")
  G7   <- e$G0[as.character(Y_mu$idr), as.character(Y_mu$idr)]
  
  list(YL = YL, G0 = e$G0, G7 = G7, X_Group = X_Group, Y_mu = Y_mu)
}

# Inverse-probability weighting for social support ## 
# ===================================================================
weight_support <- function(Z, use_lmer = FALSE, all = FALSE) {
  
  if (use_lmer) {
    # ---- GLMER branch (random intercept for cluster) ----
    covs <- if (all) {
      c("AGE", "SEX", "edu", "spouse", "DD", "Group", "Job",
        "neighbor_deg", "deg", "(1|Clus)")
    } else {
      c("AGE", "SEX", "edu", "spouse", "DD", "Job",
        "deg", "neighbor_lonely", "(1|Clus)")
    }
    fml <- as.formula(paste("Support ~", paste(covs, collapse = " + ")))
    fit <- glmer(fml, data = Z, family = "binomial",
                 control = glmerControl(optimizer = "Nelder_Mead",
                                        tolPwrss = 1e-6))
  } else {
    # ---- GLM branch (no random effects) ----
    covs <- if (all) {
      c("AGE", "SEX", "edu", "spouse", "DD", "Group", "Job",
        "neighbor_deg", "deg")
    } else {
      c("AGE", "SEX", "edu", "spouse", "DD", "Group", "Job",
        "deg", "neighbor_lonely")
    }
    fml <- as.formula(paste("Support ~", paste(covs, collapse = " + ")))
    fit <- glm(fml, data = Z, family = "binomial")
  }
  
  p     <- predict(fit, type = "response")
  mu    <- mean(Z$Support, na.rm = TRUE)
  Z$w2  <- ifelse(Z$Support == 1,
                  mu / (p + 1e-6),
                  (1 - mu) / (1 - p + 1e-6))
  Z
}

# Helper: extract emtrends from a model
# ===================================================================
.emtrend_tidy <- function(model, var) {
  summary(emtrends(model, ~Support, var = var)) %>%
    dplyr::mutate(group = dplyr::row_number()) %>%
    dplyr::select(-df) %>%
    dplyr::rename(est = !!paste0(var, ".trend")) %>%
    round(2)
}

# Separated-contacts model (lonely / non-lonely neighbour counts)
# ===================================================================
fit_separated <- function(Y, N, mnei_deg, deg_nl, mnei_lonely,
                          K = 0.5, tq = 0.005, use_lmer = FALSE) {
  
  # --- Prepare analytic sample ---
  Z <- Y %>%
    dplyr::mutate(
      deg     = deg.out,
      y       = y9,
      Support = as.integer((B.love + B.listen + B.support) > N)
    ) %>%
    dplyr::filter(!is.na(y), y7 == 0) %>%
    dplyr::select(idr, AGE, SEX, edu, spouse, D.nurse,
                  Group, Job, Support, D, y7, y, Clus) %>%
    dplyr::mutate(DD = (D + D.nurse) > 0) %>%
    dplyr::filter(!is.na(Support)) %>%
    dplyr::left_join(mnei_deg,    by = "idr") %>%
    dplyr::left_join(mnei_lonely, by = "idr") %>%
    dplyr::rename(neighbor_lonely = count) %>%
    dplyr::left_join(deg_nl,      by = "idr") %>%
    dplyr::rename(deg = count)
  
  Z$deg[is.na(Z$deg)]                         <- 0
  Z$neighbor_lonely[is.na(Z$neighbor_lonely)]  <- 0
  Z$neighbor_deg[is.na(Z$neighbor_deg)]        <- 0
  Z$deg <- (Z$deg + 1)^K - 1
  
  # --- IPW ---
  Z    <- weight_support(Z, use_lmer = use_lmer, all = FALSE)
  Z$w2 <- trimmer(Z$w2, tq)
  
  # --- Outcome model ---
  fixed <- c("AGE", "I(AGE^2)", "SEX", "edu", "spouse",
             "D", "D.nurse", "neighbor_deg", "Group", "Job",
             "neighbor_lonely * Support", "deg * Support")
  
  if (use_lmer) {
    fml <- as.formula(paste("y ~", paste(c(fixed, "(1|Clus)"),
                                         collapse = " + ")))
    m   <- glmer(fml, data = Z, weights = Z$w2, family = "binomial",
                 control = glmerControl(optimizer = "Nelder_Mead",
                                        tolPwrss = 1e-4))
  } else {
    fml <- as.formula(paste("y ~", paste(fixed, collapse = " + ")))
    m   <- glm(fml, data = Z, weights = Z$w2, family = "binomial")
  }
  
  list(
    coef       = summary(m)$coefficients,
    tab        = table(Z$Support),
    tab2       = table(Z$Support, Z$y),
    res        = .emtrend_tidy(m, "deg"),
    res.lonely = .emtrend_tidy(m, "neighbor_lonely")
  )
}

# All-contacts model (total mutual degree)
# ===================================================================
fit_allcontacts <- function(Y, N, mnei_deg,
                            tq = 0.005, use_lmer = FALSE) {
  
  Z <- Y %>%
    dplyr::mutate(
      deg     = sqrt(deg.out + 1) - 1,
      y       = y9,
      Support = as.integer((B.love + B.listen + B.support) > N)
    ) %>%
    dplyr::filter(!is.na(y), y7 == 0) %>%
    dplyr::select(idr, AGE, SEX, edu, spouse, D.nurse, B.listen, B.love,
                  Group, Job, deg, Support, D, y7, y, Clus) %>%
    dplyr::filter(!is.na(Support)) %>%
    dplyr::left_join(mnei_deg, by = "idr") %>%
    dplyr::mutate(DD = (D + D.nurse) > 0)
  
  Z$neighbor_deg[is.na(Z$neighbor_deg)] <- 0
  
  # --- IPW ---
  Z    <- weight_support(Z, use_lmer = use_lmer, all = TRUE)
  Z$w2 <- trimmer(Z$w2, tq)
  
  # --- Outcome model ---
  if (use_lmer) {
    fixed <- c("AGE", "I(AGE^2)", "SEX", "edu", "spouse", "DD",
               "Group", "Job", "neighbor_deg",
               "deg * Support", "(1|Clus)")
    fml <- as.formula(paste("y ~", paste(fixed, collapse = " + ")))
    m   <- glmer(fml, data = Z, weights = Z$w2, family = "binomial",
                 control = glmerControl(optimizer = "bobyqa"))
  } else {
    fixed <- c("AGE", "SEX", "edu", "spouse", "D",
               "Group", "Job", "neighbor_deg", "deg * Support")
    fml <- as.formula(paste("y ~", paste(fixed, collapse = " + ")))
    m   <- glm(fml, data = Z, weights = Z$w2, family = "binomial")
  }
  
  list(
    coef = summary(m)$coefficients,
    tab  = table(Z$Support),
    tab2 = table(Z$Support, Z$y),
    res  = .emtrend_tidy(m, "deg")
  )
}

# Collect marginal betas across thresholds
# ===================================================================
extract_betas <- function(sep, allc = NULL) {
  B <- list(
    ls_d = do.call(rbind, lapply(sep, \(x) x$res[1, ])),
    su_d = do.call(rbind, lapply(sep, \(x) x$res[2, ])),
    ls_C = do.call(rbind, lapply(sep, \(x) x$res.lonely[1, ])),
    su_C = do.call(rbind, lapply(sep, \(x) x$res.lonely[2, ]))
  )
  if (!is.null(allc)) {
    B$deg_l <- do.call(rbind, lapply(allc, \(x) x$res[1, ]))
    B$deg_m <- do.call(rbind, lapply(allc, \(x) x$res[2, ]))
  }
  B
}
