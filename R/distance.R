# ---------- 1) Distances and maxima (helpers) ----------
#' @noRd
kendall_distance <- function(p, q) {
  K <- length(p)
  idx <- utils::combn(K, 2)
  sum( (p[idx[1,]] - p[idx[2,]]) * (q[idx[1,]] - q[idx[2,]]) < 0 )
}

#' @noRd
footrule_distance <- function(p, q) sum(abs(p - q))

#' @noRd
rho_distance <- function(p, q) sum((p - q)^2)

#' @noRd
dmax_for <- function(K, method = c("kendall", "footrule", "rho")) {
  method <- match.arg(method)
  switch(method,
         "kendall"  = K*(K-1)/2,
         "footrule" = if (K %% 2 == 0) K^2/2 else (K^2 - 1)/2,
         "rho"      = K*(K^2 - 1)/3)
}

#' @noRd
dist_fun_for <- function(method = c("kendall","footrule","rho")) {
  method <- match.arg(method)
  switch(method,
         "kendall"  = kendall_distance,
         "footrule" = footrule_distance,
         "rho"      = rho_distance)
}


# ---------- 2) Main function ----------
#' Estimate a Cross-Group Rank Gap (U-statistic)
#'
#' Computes the average pairwise rank-distance between two groups and normalizes
#' it by the method-specific maximum, yielding \eqn{\theta \in [0,1]}. Supports
#' analytic or bootstrap standard errors.
#'
#' @param df A data frame containing rankings of the same \eqn{K} items.
#' @param group_col Name of the column with the grouping variable (must have exactly two groups).
#' @param rank_cols Character vector of column names holding the item ranks (1 = best).
#' @param method Distance to use: \code{"kendall"}, \code{"footrule"}, or \code{"rho"}.
#' @param se SE method: \code{"bootstrap"} or \code{"analytic"}.
#' @param B Number of bootstrap resamples if \code{se = "bootstrap"}.
#' @param seed Optional random seed for reproducibility.
#'
#' @return An object of class \code{rank_gap} with fields:
#' \itemize{
#' \item \code{theta_hat}: normalized gap in \eqn{[0,1]}.
#' \item \code{se}: standard error.
#' \item \code{ci95}: 95\% interval.
#' \item \code{nT}, \code{nC}, \code{K}, \code{d_max}, \code{method}, \code{type}, \code{estimator}.
#' }
#'
#' @importFrom dplyr filter select group_by summarise across bind_rows all_of %>%
#' @importFrom rlang .data
#' @importFrom utils combn
#' @importFrom stats var quantile sd
#'
#' @examples
#' set.seed(1)
#' K <- 5; items <- c("A","B","C","D","E")
#' # controls: random permutations
#' Cmat <- as.data.frame(t(replicate(80, sample(1:K, K))))
#' names(Cmat) <- items
#' # treated: like controls but bias item C to be worse
#' Tmat <- as.data.frame(t(replicate(80, {
#'   p <- sample(1:K, K)
#'   p[which(items=="C")] <- pmin(K, p[which(items=="C")] + sample(0:2, 1))
#'   rank(-rank(p, ties.method="first"))
#' })))
#' names(Tmat) <- items
#' # build data (no tibble needed)
#' df <- dplyr::bind_rows(
#'   data.frame(group="T", Tmat, check.names = FALSE),
#'   data.frame(group="C", Cmat, check.names = FALSE)
#' )
#' res <- est_rank_gap(df, group_col="group", rank_cols=items,
#'                     method="kendall", se="bootstrap", B=200, seed=42)
#' print(res)

#' @export
est_rank_gap <- function(df, group_col, rank_cols,
                         method = c("kendall","footrule","rho"),
                         se = c("bootstrap","analytic"),
                         B = 1000, seed = NULL) {

  method <- match.arg(method)
  se     <- match.arg(se)

  g <- df[[group_col]]
  stopifnot(all(!is.na(g)))
  levs <- unique(g)
  stopifnot(length(levs) == 2)
  gT <- levs[1]; gC <- levs[2]  # will keep this order

  X <- as.matrix(df %>% dplyr::filter(.data[[group_col]] == gT) %>%
                   dplyr::select(dplyr::all_of(rank_cols)))
  Y <- as.matrix(df %>% dplyr::filter(.data[[group_col]] == gC) %>%
                   dplyr::select(dplyr::all_of(rank_cols)))

  nT <- nrow(X); nC <- nrow(Y); K <- ncol(X)
  stopifnot(K == ncol(Y))

  dfun  <- dist_fun_for(method)
  dmax  <- dmax_for(K, method)

  # distance matrix (nT x nC)
  D <- matrix(0, nT, nC)
  for (i in 1:nT) for (j in 1:nC) D[i, j] <- dfun(X[i, ], Y[j, ])

  U_hat_raw <- mean(D)
  theta_hat <- U_hat_raw / dmax

  if (se == "analytic") {
    row_means <- rowMeans(D)
    col_means <- colMeans(D)
    var_hat   <- stats::var(row_means)/nT + stats::var(col_means)/nC
    se_hat    <- sqrt(var_hat) / dmax
    ci        <- theta_hat + c(-1,1)*1.96*se_hat

    out <- list(method = method, type = "cross_pairs",
                nT = nT, nC = nC, K = K, d_max = dmax,
                theta_hat = theta_hat, se = se_hat, ci95 = ci,
                estimator = "U-statistic (analytic SE)")
  } else {
    if (!is.null(seed)) set.seed(seed)
    boots <- numeric(B)
    for (b in 1:B) {
      iT <- sample.int(nT, nT, TRUE)
      iC <- sample.int(nC, nC, TRUE)
      Db <- D[iT, iC, drop = FALSE]
      boots[b] <- mean(Db) / dmax
    }
    se_hat <- sd(boots)
    ci     <- stats::quantile(boots, c(0.025, 0.975), names = FALSE)

    out <- list(method = method, type = "cross_pairs",
                nT = nT, nC = nC, K = K, d_max = dmax,
                theta_hat = theta_hat, se = se_hat, ci95 = ci,
                estimator = paste0("U-statistic (bootstrap B=", B, ")"),
                boot_samples = boots)
  }
  class(out) <- "rank_gap"
  out
}


# ---------- 3) S3 print method ----------
#' Print a rank_gap object
#' @param x A \code{rank_gap} object.
#' @param ... Passed to methods.
#' @export
#' @method print rank_gap
print.rank_gap <- function(x, ...) {
  cat(sprintf("Method: %s | Type: %s | K=%d | nT=%d, nC=%d\n",
              x$method, x$type, x$K, x$nT, x$nC))
  cat(sprintf("d_max = %.2f\n", x$d_max))
  cat(sprintf("theta_hat (normalized 0-1) = %.4f\n", x$theta_hat))
  cat(sprintf("SE = %.4f | 95%% CI = [%.4f, %.4f]\n",
              x$se, x$ci95[1], x$ci95[2]))
}


# ---------- 4) Optional consensus variant (export or not as you wish) ----------
#' @noRd
est_rank_gap_consensus <- function(df, group_col, rank_cols,
                                   method = c("footrule","rho"),
                                   B = 1000, seed = NULL) {
  method <- match.arg(method)
  dfun   <- dist_fun_for(method)
  K      <- length(rank_cols)
  dmax   <- dmax_for(K, method)

  avg_rank <- df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(rank_cols), mean), .groups = "drop")
  stopifnot(nrow(avg_rank) == 2)

  rT <- as.numeric(avg_rank[1, rank_cols])
  rC <- as.numeric(avg_rank[2, rank_cols])
  theta_hat <- dfun(rT, rC) / dmax

  g <- df[[group_col]]
  levs <- unique(g); gT <- levs[1]; gC <- levs[2]
  X <- df %>% dplyr::filter(.data[[group_col]] == gT) %>% dplyr::select(dplyr::all_of(rank_cols)) %>% as.matrix()
  Y <- df %>% dplyr::filter(.data[[group_col]] == gC) %>% dplyr::select(dplyr::all_of(rank_cols)) %>% as.matrix()
  nT <- nrow(X); nC <- nrow(Y)

  if (!is.null(seed)) set.seed(seed)
  boots <- numeric(B)
  for (b in 1:B) {
    rT_b <- colMeans(X[sample.int(nT, nT, TRUE), , drop = FALSE])
    rC_b <- colMeans(Y[sample.int(nC, nC, TRUE), , drop = FALSE])
    boots[b] <- dfun(rT_b, rC_b) / dmax
  }
  se_hat <- sd(boots)
  ci     <- stats::quantile(boots, c(0.025, 0.975), names = FALSE)

  out <- list(method = method, type = "consensus_borda",
              nT = nT, nC = nC, K = K, d_max = dmax,
              theta_hat = theta_hat, se = se_hat, ci95 = ci,
              estimator = paste0("Consensus (bootstrap B=", B, ")"))
  class(out) <- "rank_gap"
  out
}
