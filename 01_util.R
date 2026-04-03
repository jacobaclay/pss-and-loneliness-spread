library(igraph)
library(randnet)
library(dplyr)

# ---------------------------------------------------------------------------
# Network statistics
# ---------------------------------------------------------------------------
get_network_stat <- function(A, stat = "degree") {
  ids <- as.numeric(colnames(A))
  
  if (stat == "degree") {
    data.frame(idr = ids, deg.in = colSums(A), deg.out = rowSums(A))
    
  } else if (stat == "constraint") {

    G <- graph_from_adjacency_matrix(A)
    data.frame(idr = ids, Constraint = igraph::constraint(G))
 
  } else if (stat == "transitivity") {
    G <- graph_from_adjacency_matrix(A)
    data.frame(idr = ids,
               transitivity = igraph::transitivity(G, type = "local"))
  }
}

# ---------------------------------------------------------------------------
# Expand adjacency matrix to include new nodes (with zero ties)
# ---------------------------------------------------------------------------
expand_matrix <- function(U, A) {
  new <- setdiff(as.character(U$idr), colnames(A))
  if (length(new) == 0) return(A)
  n <- nrow(A) + length(new)
  A2 <- matrix(0, n, n)
  nms <- c(rownames(A), new)
  rownames(A2) <- colnames(A2) <- nms
  A2[rownames(A), colnames(A)] <- A
  A2
}

# ---------------------------------------------------------------------------
# Winsorise extreme weights
# ---------------------------------------------------------------------------
trimmer <- function(x, q = 0.0125) {
  lo <- quantile(x, q, na.rm = TRUE)
  hi <- quantile(x, 1 - q, na.rm = TRUE)
  pmax(pmin(x, hi), lo)
}

# ---------------------------------------------------------------------------
# Mutual ties to lonely (y7==1) or non-lonely (y7==0) neighbours
# ---------------------------------------------------------------------------
find_mutual_ties <- function(Y, Net, lonely = TRUE) {
  Yf <- Y %>%
    dplyr::select(idr, y7) %>%
    dplyr::filter(idr %in% as.numeric(colnames(Net)), !is.na(y7))
  
  target <- Yf %>%
    dplyr::filter(y7 == as.integer(lonely)) %>%
    dplyr::pull(idr)
  
  ids <- intersect(rownames(Net), colnames(Net))
  Nm  <- (Net[ids, ids] > 0) & t(Net[ids, ids] > 0)
  sub <- (Nm * 1)[as.character(Yf$idr), as.character(target), drop = FALSE]
  
  data.frame(idr = as.numeric(rownames(sub)), count = rowSums(sub))
}

# ---------------------------------------------------------------------------
# Average neighbor degree (power-transformed)
# ---------------------------------------------------------------------------
find_avg_neighbor_degree <- function(Y, Net, K = 0.5) {
  Yf <- Y %>%
    dplyr::mutate(deg = deg.out) %>%
    dplyr::select(idr, y7, y9, deg) %>%
    dplyr::filter(idr %in% as.numeric(colnames(Net)), !is.na(y7))
  
  A_sym <- ((Net + t(Net)) > 0)
  A_sub <- A_sym[as.character(Yf$idr), as.character(Yf$idr)]
  
  avg <- as.numeric(((A_sub %*% Yf$deg + 1)^K - 1) / rowSums(A_sub))
  avg[!is.finite(avg)] <- NA
  
  data.frame(idr = Yf$idr, neighbor_deg = avg)
}
