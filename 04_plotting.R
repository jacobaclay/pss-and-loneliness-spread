source('config.R')

# ---------------------------------------------------------------------------
# Step-plot helpers
# ---------------------------------------------------------------------------
.step <- function(i, d, col, lwd = 3) {
  
  dash <- ifelse(d$sig[i + 1] == 0, 3, 1)
  lines(i:(i + 1), rep(d$est[i + 1], 2),
        type = "s", lwd = lwd, lty = dash, col = col)
}

.step_ci <- function(i, d, col, lwd = 4) {
  lines(i:(i + 1), rep(d$est[i + 1], 2) - 0.04, lwd = 1.2, col = col)
  lines(i:(i + 1), rep(d$est[i + 1], 2) + 0.04, lwd = 1.2, col = col)
  .step(i, d, col, lwd)
}

add_sig <- function(df, z = 1.9) {
  df %>% dplyr::mutate(sig = as.integer(SE * z < abs(est)))
}

# ---------------------------------------------------------------------------
# Main marginal-beta figure
# ---------------------------------------------------------------------------
plot_marginal <- function(B, ylim = c(-2, 2),
                          labels = paste0("S", seq_along(CONFIG$thresholds))) {
  d  <- lapply(B[c("ls_d", "su_d", "ls_C", "su_C")], add_sig)
  d2 <- lapply(B[c("deg_l", "deg_m")], add_sig, z = 1.8)
  nr <- nrow(d[[1]]); nr1 <- nr - 1
  
  plot(0:nr1, rep(0, nr), ylim = ylim, col = "white",
       xlab = "", ylab = expression(paste("Marginal ", beta)),
       xaxt = "n", xlim = c(0, nr))
  
  lapply(0:nr1, \(i) .step(i, d$su_d, rgb(0, 0, .75, .35), 4))
  lapply(0:nr1, \(i) .step(i, d$su_C, rgb(.9, 0, 0, .35),  4))
  lapply(0:nr1, \(i) .step(i, d$ls_d, "blue"))
  lapply(0:nr1, \(i) .step(i, d$ls_C, "red"))
  lapply(0:nr1, \(i) .step_ci(i, d2$deg_m, "grey70", 3))
  lapply(0:nr1, \(i) .step_ci(i, d2$deg_l, "black",  2))
  
  lines(0:nr, rep(0, nr + 1), col = rgb(0, 0, 0, .1), lwd = 9)
  for (v in seq_len(nr1)) {
    abline(v = v, lwd = 8, col = "grey25")
    abline(v = v, lwd = 2, col = "white")
  }
  axis(1, 0.5 + 0:nr1, labels = labels)
}

# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------
plot_legend <- function() {
  par(mfrow = c(1, 1))
  plot(0:5, 0:5, col = "white", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "")
  legend(1, 3,
         legend = c("Non-lonely Contacts", "Lonely Contacts"),
         title  = "Supported (PSS+)", title.font = 2,
         lwd = 5, col = c(rgb(0, 0, .5, .25), rgb(.5, 0, 0, .25)),
         bty = "n")
  legend(2.5, 3,
         legend = c("Non-lonely Contacts", "Lonely Contacts"),
         title  = "Lacks Support (PSS-)", title.font = 2,
         lwd = 5, col = c("blue", "red"), bty = "n")
  legend(1.25, 2,
         legend = c("Significant at 95%", "Non-significant"),
         lwd = c(5, 3), lty = c(1, 3), text.font = 3,
         bty = "n", horiz = TRUE)
}

# ---------------------------------------------------------------------------
# Diagnostic: degree transformations
# ---------------------------------------------------------------------------
check_sqrt <- function(deg) {
  par(mfrow = c(1, 3))
  hist(deg, breaks = 15, main = "Untransformed",          cex.main = 1.2)
  hist(sqrt(deg), breaks = 15, main = "sqrt(x)",
       col = rgb(.5, .5, 1, .2), cex.main = 1.2)
  hist(sqrt(deg + 1) - 1, breaks = 12, main = "sqrt(x+1)-1",
       col = rgb(0, 0, 1, .2), cex.main = 1.2)
}