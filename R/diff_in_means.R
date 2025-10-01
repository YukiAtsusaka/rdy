#' Estimate ranking effects: weighted, unit-normalized [-1,1], and standardized
#'
#' For each item in \code{items}, constructs:
#' 1) Average rank of the target (1 = best), reversed and scaled to unit range:
#'    \eqn{\Delta_{\text{unit}} = -\Delta_{\text{rank}}/(J-1) \in [-1,1]}.
#' 2) Pairwise wins vs. others (binary).
#' 3) Top-k indicators for the target (k = 1..J-1) (binary).
#' 4) Marginal rank indicators for the target (r = 1..J) (binary).
#'
#' Estimation is difference-in-means via \code{estimatr::lm_robust} (optionally weighted).
#' For each effect, returns:
#'   - Raw ATE (\code{*_raw})
#'   - \strong{Unit-normalized} ATE in [-1,1] (\code{*_unit})
#'   - \strong{Standardized} ATE dividing \code{*_unit} by the control SD on the same scale (\code{*_std})
#'
#' @param data A data.frame/tibble containing \code{items}, \code{treat}, and optionally \code{weights}.
#' @param items Character vector of ranking columns (1 = top/best). Length \(\ge 2\).
#' @param treat String; binary treatment column (0/1, logical, or factor that can be mapped to 0/1).
#' @param weights NULL, a numeric vector of length \code{nrow(data)}, or a single string naming a weight column.
#' @param se_type SE type passed to \code{estimatr::lm_robust} (default "HC2").
#' @param clusters Optional string naming a clustering variable in \code{data} (forwarded to \code{lm_robust}).
#' @param na_action "drop" (default) to listwise drop on used columns, or "none".
#' @param top_k_max Optional integer; cap for top-k block (default \code{length(items)-1}).
#' @param return_models Logical; if TRUE, also returns (placeholder) model list. Default FALSE.
#' @param ... Extra args passed to \code{estimatr::lm_robust}.
#'
#' @return Tibble with rows for all effects across all items. Columns include:
#' \itemize{
#'  \item \code{estimate_raw, std.error_raw, conf.low_raw, conf.high_raw}
#'  \item \code{estimate_unit, std.error_unit, conf.low_unit, conf.high_unit}  (in [-1,1])
#'  \item \code{estimate_std, std.error_std, conf.low_std, conf.high_std}      (standardized by control SD)
#'  \item \code{sd_control_unit} (the control-group SD on the unit scale used to standardize)
#'  \item \code{p.value}, \code{outcome}, \code{effect_type}, \code{item}, \code{J}
#'  \item diagnostics: \code{n, n_treat, n_control, w_sum, w_treat_sum, w_control_sum}
#' }
#' If \code{return_models=TRUE}, returns a list with \code{results} and a placeholder \code{models}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dt <- data.frame(
#'   treated = rbinom(300, 1, 0.5),
#'   A = sample(1:4, 300, TRUE),
#'   B = sample(1:4, 300, TRUE),
#'   C = sample(1:4, 300, TRUE),
#'   w = runif(300, 0.5, 1.5)
#' )
#' res <- diff_in_means(dt, items = c("A","B","C"), treat = "treated", weights = "w")
#' }
#' @export
diff_in_means <- function(data,
                       items,
                       treat,
                       weights = NULL,
                       se_type = "HC2",
                       clusters = NULL,
                       na_action = c("drop", "none"),
                       top_k_max = NULL,
                       return_models = FALSE,
                       ...) {
  # --- namespaces ---
  requireNamespace("estimatr", quietly = TRUE)
  requireNamespace("broom",    quietly = TRUE)
  requireNamespace("dplyr",    quietly = TRUE)
  requireNamespace("purrr",    quietly = TRUE)
  requireNamespace("rlang",    quietly = TRUE)
  requireNamespace("tibble",   quietly = TRUE)

  na_action <- match.arg(na_action)

  # --- checks ---
  stopifnot(is.data.frame(data))
  if (!is.character(items) || length(items) < 2L)
    stop("`items` must be a character vector with at least two ranking columns.")
  if (!is.character(treat) || length(treat) != 1L)
    stop("`treat` must be a single column name (string).")
  if (!all(items %in% names(data))) stop("Some `items` not found in `data`.")
  if (!treat %in% names(data))      stop("`treat` not found in `data`.")
  if (!is.null(clusters) && !clusters %in% names(data))
    stop("`clusters` not found in `data`.")

  # --- weights: vector or column ---
  w <- NULL
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1L) {
      if (!weights %in% names(data)) stop("`weights` column not found in `data`.")
      w <- data[[weights]]
    } else if (is.numeric(weights)) {
      if (length(weights) != nrow(data))
        stop("If `weights` is numeric, it must have length nrow(data).")
      w <- weights
    } else {
      stop("`weights` must be NULL, a single column name, or a numeric vector.")
    }
  }

  used_cols <- unique(c(items, treat, clusters, if (is.character(weights)) weights))
  df <- data[, used_cols, drop = FALSE]

  if (na_action == "drop") {
    keep <- stats::complete.cases(df)
    df <- df[keep, , drop = FALSE]
    if (!is.null(w)) w <- w[keep]
  }

  # --- treatment to {0,1} ---
  D <- df[[treat]]
  if (is.logical(D)) D <- as.integer(D)
  if (is.factor(D))  D <- as.integer(D) - 1L
  if (!all(D %in% c(0, 1)))
    stop("`treat` must be binary/codable to {0,1} (0 = control, 1 = treated).")

  # diagnostics
  n       <- nrow(df)
  n_treat <- sum(D == 1, na.rm = TRUE)
  n_ctrl  <- sum(D == 0, na.rm = TRUE)
  if (!is.null(w)) {
    w_treat_sum <- sum(w[D == 1], na.rm = TRUE)
    w_ctrl_sum  <- sum(w[D == 0], na.rm = TRUE)
    w_sum       <- sum(w, na.rm = TRUE)
  } else {
    w_treat_sum <- n_treat
    w_ctrl_sum  <- n_ctrl
    w_sum       <- n
  }

  J <- length(items)
  if (is.null(top_k_max)) top_k_max <- J - 1L
  top_k_max <- max(1L, min(as.integer(top_k_max), J - 1L))

  # --- helpers ---
  w_or_ones <- function(x) if (is.null(x)) rep(1, nrow(df)) else x

  wmean <- function(y, ww) {
    ww <- w_or_ones(ww)
    sum(ww * y, na.rm = TRUE) / sum(ww, na.rm = TRUE)
  }
  wsd_pop <- function(y, ww) {
    ww <- w_or_ones(ww)
    mu <- wmean(y, ww)
    v  <- sum(ww * (y - mu)^2, na.rm = TRUE) / sum(ww, na.rm = TRUE)
    sqrt(pmax(v, 0))
  }

  # run weighted diff-in-means via lm_robust on vector y
  run_dim <- function(y) {
    fit <- estimatr::lm_robust(
      y ~ D,
      data    = df,
      weights = w,
      se_type = se_type,
      clusters = if (!is.null(clusters)) df[[clusters]] else NULL,
      ...
    )
    broom::tidy(fit, conf.int = TRUE)
  }

  # linear transform of tidy (x -> a + b x), ordering CI properly for b<0
  lin_xform <- function(tb, a, b) {
    lo <- a + b * tb$conf.low
    hi <- a + b * tb$conf.high
    tb$estimate  <- a + b * tb$estimate
    tb$std.error <- abs(b) * tb$std.error
    tb$conf.low  <- pmin(lo, hi)
    tb$conf.high <- pmax(lo, hi)
    tb
  }

  # assemble a block given a target
  by_target <- function(target_item) {
    other_items <- setdiff(items, target_item)
    rank_t     <- df[[target_item]]

    # ---- Average rank block ----
    # Raw Δ_rank in [-(J-1), (J-1)]
    avg_raw <- run_dim(rank_t) |>
      dplyr::filter(.data$term == "D") |>
      dplyr::mutate(
        outcome     = paste0("Avg: ", target_item),
        effect_type = "average rank",
        item        = target_item,
        J           = J
      )

    # Unit-normalized ([-1,1]): Δ_unit = -Δ_rank/(J-1)
    b_avg <- -1 / (J - 1)
    avg_unit <- lin_xform(dplyr::rename(avg_raw,
                                        estimate_raw  = estimate,
                                        std.error_raw = std.error,
                                        conf.low_raw  = conf.low,
                                        conf.high_raw = conf.high),
                          a = 0, b = b_avg)

    # control SD on the same unit scale: sd(rank|D=0)/(J-1)
    ctrl_mask <- (D == 0)
    sd_ctrl_avg_unit <- wsd_pop(rank_t[ctrl_mask], if (is.null(w)) NULL else w[ctrl_mask]) / (J - 1)
    if (is.na(sd_ctrl_avg_unit) || sd_ctrl_avg_unit <= .Machine$double.eps) {
      avg_unit <- avg_unit |>
        dplyr::mutate(
          estimate_std  = NA_real_,
          std.error_std = NA_real_,
          conf.low_std  = NA_real_,
          conf.high_std = NA_real_,
          sd_control_unit = sd_ctrl_avg_unit
        )
    } else {
      avg_unit <- avg_unit |>
        dplyr::mutate(
          estimate_std  = estimate / sd_ctrl_avg_unit,
          std.error_std = std.error / sd_ctrl_avg_unit,
          conf.low_std  = conf.low / sd_ctrl_avg_unit,
          conf.high_std = conf.high / sd_ctrl_avg_unit,
          sd_control_unit = sd_ctrl_avg_unit
        )
    }

    # keep tidy columns
    res_avg <- avg_unit |>
      dplyr::select(
        estimate_raw, std.error_raw, conf.low_raw, conf.high_raw,
        estimate, std.error, conf.low, conf.high,
        estimate_std, std.error_std, conf.low_std, conf.high_std,
        sd_control_unit, p.value,
        outcome, effect_type, item, J
      ) |>
      dplyr::rename(
        estimate_unit  = estimate,
        std.error_unit = std.error,
        conf.low_unit  = conf.low,
        conf.high_unit = conf.high
      )

    # ---- Pairwise wins (binary) ----
    pair_list <- purrr::imap(other_items, function(oi, j) {
      y <- as.integer(rank_t < df[[oi]])
      tb <- run_dim(y) |>
        dplyr::filter(.data$term == "D") |>
        dplyr::mutate(
          outcome     = paste0("v. ", oi),
          effect_type = "pairwise ranking",
          item        = target_item,
          J           = J
        )
      # raw Δ in [-1,1]; unit = raw (b=1)
      tb2 <- dplyr::rename(tb,
                           estimate_raw  = estimate,
                           std.error_raw = std.error,
                           conf.low_raw  = conf.low,
                           conf.high_raw = conf.high)
      tb2 <- lin_xform(tb2, a = 0, b = 1)

      y0  <- y[ctrl_mask]
      sd0 <- wsd_pop(y0, if (is.null(w)) NULL else w[ctrl_mask])
      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        tb2 <- tb2 |>
          dplyr::mutate(
            estimate_std  = NA_real_,
            std.error_std = NA_real_,
            conf.low_std  = NA_real_,
            conf.high_std = NA_real_,
            sd_control_unit = sd0
          )
      } else {
        tb2 <- tb2 |>
          dplyr::mutate(
            estimate_std  = estimate / sd0,
            std.error_std = std.error / sd0,
            conf.low_std  = conf.low / sd0,
            conf.high_std = conf.high / sd0,
            sd_control_unit = sd0
          )
      }

      tb2 |>
        dplyr::select(
          estimate_raw, std.error_raw, conf.low_raw, conf.high_raw,
          estimate, std.error, conf.low, conf.high,
          estimate_std, std.error_std, conf.low_std, conf.high_std,
          sd_control_unit, p.value,
          outcome, effect_type, item, J
        ) |>
        dplyr::rename(
          estimate_unit  = estimate,
          std.error_unit = std.error,
          conf.low_unit  = conf.low,
          conf.high_unit = conf.high
        )
    })
    res_pair <- dplyr::bind_rows(pair_list)

    # ---- Top-k (binary) ----
    top_list <- purrr::imap(seq_len(top_k_max), function(k, j) {
      y <- as.integer(rank_t <= k)
      tb <- run_dim(y) |>
        dplyr::filter(.data$term == "D") |>
        dplyr::mutate(
          outcome     = paste0("Top-", k),
          effect_type = "top-k ranking",
          item        = target_item,
          J           = J
        )
      tb2 <- dplyr::rename(tb,
                           estimate_raw  = estimate,
                           std.error_raw = std.error,
                           conf.low_raw  = conf.low,
                           conf.high_raw = conf.high)
      tb2 <- lin_xform(tb2, a = 0, b = 1)

      y0  <- y[ctrl_mask]
      sd0 <- wsd_pop(y0, if (is.null(w)) NULL else w[ctrl_mask])
      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        tb2 <- tb2 |>
          dplyr::mutate(
            estimate_std  = NA_real_,
            std.error_std = NA_real_,
            conf.low_std  = NA_real_,
            conf.high_std = NA_real_,
            sd_control_unit = sd0
          )
      } else {
        tb2 <- tb2 |>
          dplyr::mutate(
            estimate_std  = estimate / sd0,
            std.error_std = std.error / sd0,
            conf.low_std  = conf.low / sd0,
            conf.high_std = conf.high / sd0,
            sd_control_unit = sd0
          )
      }

      tb2 |>
        dplyr::select(
          estimate_raw, std.error_raw, conf.low_raw, conf.high_raw,
          estimate, std.error, conf.low, conf.high,
          estimate_std, std.error_std, conf.low_std, conf.high_std,
          sd_control_unit, p.value,
          outcome, effect_type, item, J
        ) |>
        dplyr::rename(
          estimate_unit  = estimate,
          std.error_unit = std.error,
          conf.low_unit  = conf.low,
          conf.high_unit = conf.high
        )
    })
    res_top <- dplyr::bind_rows(top_list)

    # ---- Marginal rank (binary) ----
    marg_list <- purrr::imap(seq_len(J), function(r, j) {
      y <- as.integer(rank_t == r)
      tb <- run_dim(y) |>
        dplyr::filter(.data$term == "D") |>
        dplyr::mutate(
          outcome     = paste0("Ranked ", r),
          effect_type = "marginal ranking",
          item        = target_item,
          J           = J
        )
      tb2 <- dplyr::rename(tb,
                           estimate_raw  = estimate,
                           std.error_raw = std.error,
                           conf.low_raw  = conf.low,
                           conf.high_raw = conf.high)
      tb2 <- lin_xform(tb2, a = 0, b = 1)

      y0  <- y[ctrl_mask]
      sd0 <- wsd_pop(y0, if (is.null(w)) NULL else w[ctrl_mask])
      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        tb2 <- tb2 |>
          dplyr::mutate(
            estimate_std  = NA_real_,
            std.error_std = NA_real_,
            conf.low_std  = NA_real_,
            conf.high_std = NA_real_,
            sd_control_unit = sd0
          )
      } else {
        tb2 <- tb2 |>
          dplyr::mutate(
            estimate_std  = estimate / sd0,
            std.error_std = std.error / sd0,
            conf.low_std  = conf.low / sd0,
            conf.high_std = conf.high / sd0,
            sd_control_unit = sd0
          )
      }

      tb2 |>
        dplyr::select(
          estimate_raw, std.error_raw, conf.low_raw, conf.high_raw,
          estimate, std.error, conf.low, conf.high,
          estimate_std, std.error_std, conf.low_std, conf.high_std,
          sd_control_unit, p.value,
          outcome, effect_type, item, J
        ) |>
        dplyr::rename(
          estimate_unit  = estimate,
          std.error_unit = std.error,
          conf.low_unit  = conf.low,
          conf.high_unit = conf.high
        )
    })
    res_marg <- dplyr::bind_rows(marg_list)

    dplyr::bind_rows(res_avg, res_pair, res_top, res_marg)
  }

  # run for all targets
  results <- purrr::map_dfr(items, by_target)

  # add global diagnostics
  results <- results |>
    dplyr::mutate(
      n             = n,
      n_treat       = n_treat,
      n_control     = n_ctrl,
      w_sum         = w_sum,
      w_treat_sum   = w_treat_sum,
      w_control_sum = w_ctrl_sum
    ) |>
    tibble::as_tibble()

  if (!return_models) return(results)

  # placeholder to keep API parity
  list(results = results, models = NULL)
}
