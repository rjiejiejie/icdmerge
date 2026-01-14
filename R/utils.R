# internal utilities ------------------------------------------------------

row_zscore <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  row_means <- rowMeans(X)
  Xc <- X - row_means
  row_sds <- sqrt(rowMeans(Xc^2))
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  Xz <- Xc / row_sds
  Xz[is.na(Xz)] <- 0
  Xz
}

cocluster_mat <- function(cluster, k) {
  n <- length(cluster)
  Z <- Matrix::sparseMatrix(i = 1:n, j = cluster, x = 1, dims = c(n, k))
  C <- Matrix::tcrossprod(Z)
  C <- as.matrix(C)
  diag(C) <- 0
  C
}

calc_pac <- function(cons_prob, x1 = 0.1, x2 = 0.9) {
  vals <- cons_prob[upper.tri(cons_prob, diag = FALSE)]
  if (length(vals) == 0) return(NA_real_)
  mean(vals > x1 & vals < x2)
}

sanitize_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

calc_homo_q <- function(cons_prob, q = 0.1) {
  vals <- cons_prob[upper.tri(cons_prob, diag = FALSE)]
  if (length(vals) == 0) return(NA_real_)
  as.numeric(stats::quantile(vals, probs = q, na.rm = TRUE, names = FALSE, type = 7))
}
