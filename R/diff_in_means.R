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
#'   - \strong{Standardized} ATE dividing \code{_unit} by the control SD on the same scale (\code{*_std})
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
#' @return Tibble with rows for all effects across all items.
#' @importFrom rlang .data
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
  requireNamespace("estimatr", quietly = TRUE)
  requireNamespace("broom",    quietly = TRUE)
  requireNamespace("dplyr",    quietly = TRUE)
  requireNamespace("purrr",    quietly = TRUE)
  requireNamespace("tibble",   quietly = TRUE)
  requireNamespace("tidyr",    quietly = TRUE)

  na_action <- match.arg(na_action)
  stopifnot(is.data.frame(data))

  if (!is.character(items) || length(items) < 2L)
    stop("`items` must be a character vector with at least two ranking columns.")
  if (!is.character(treat) || length(treat) != 1L)
    stop("`treat` must be a single column name (string).")
  if (!all(items %in% names(data))) stop("Some `items` not found in `data`.")
  if (!treat %in% names(data))      stop("`treat` not found in `data`.")
  if (!is.null(clusters) && !clusters %in% names(data))
    stop("`clusters` not found in `data`.")

  # Weights
  w <- NULL
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1L) {
      if (!weights %in% names(data)) stop("`weights` column not found in `data`.")
      w <- data[[weights]]
    } else if (is.numeric(weights)) {
      if (length(weights) != nrow(data))
        stop("If `weights` is numeric, it must have length nrow(data).")
      w <- weights
    } else stop("`weights` must be NULL, a single column name, or a numeric vector.")
  }

  used_cols <- unique(c(items, treat, clusters, if (is.character(weights)) weights))
  df <- data[, used_cols, drop = FALSE]
  if (na_action == "drop") {
    keep <- stats::complete.cases(df)
    df <- df[keep, , drop = FALSE]
    if (!is.null(w)) w <- w[keep]
  }

  # Treatment coding
  D <- df[[treat]]
  if (is.logical(D)) D <- as.integer(D)
  if (is.factor(D))  D <- as.integer(D) - 1L
  if (!all(D %in% c(0, 1)))
    stop("`treat` must be binary/codable to {0,1} (0 = control, 1 = treated).")

  # Diagnostics
  n       <- nrow(df)
  n_treat <- sum(D == 1, na.rm = TRUE)
  n_ctrl  <- sum(D == 0, na.rm = TRUE)

  if (is.null(w)) {
    w_vec <- rep(1, n)
    w_treat_sum <- n_treat
    w_ctrl_sum  <- n_ctrl
    w_sum       <- n
  } else {
    w_vec <- w
    w_treat_sum <- sum(w_vec[D == 1], na.rm = TRUE)
    w_ctrl_sum  <- sum(w_vec[D == 0], na.rm = TRUE)
    w_sum       <- sum(w_vec, na.rm = TRUE)
  }

  # Settings
  J <- length(items)
  if (is.null(top_k_max)) top_k_max <- J - 1L
  top_k_max <- max(1L, min(as.integer(top_k_max), J - 1L))
  ctrl_mask <- (D == 0)

  # Utilities
  wmean_masked <- function(y, ww, mask) {
    if (!any(mask, na.rm = TRUE)) return(NA_real_)
    num <- sum(ww[mask] * y[mask], na.rm = TRUE)
    den <- sum(ww[mask],           na.rm = TRUE)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    num / den
  }

  wsd_masked_fsc <- function(y, ww, mask) {
    idx <- which(mask & !is.na(y) & !is.na(ww))
    if (length(idx) == 0) return(NA_real_)
    w  <- ww[idx]
    y2 <- y[idx]
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) return(NA_real_)
    w  <- w / sw
    mu <- sum(w * y2)
    denom <- 1 - sum(w^2)
    if (denom <= .Machine$double.eps) return(NA_real_)
    v <- sum(w * (y2 - mu)^2) / denom
    sqrt(pmax(v, 0))
  }

  fit_dim <- function(y) {
    estimatr::lm_robust(
      y ~ D, data = df,
      weights  = if (!is.null(weights)) w_vec else NULL,
      se_type  = se_type,
      clusters = if (!is.null(clusters)) df[[clusters]] else NULL,
      ...
    )
  }

  to_rank_num <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    if (is.factor(x))  return(suppressWarnings(as.numeric(as.character(x))))
    if (is.character(x)) return(suppressWarnings(as.numeric(x)))
    stop("Rank columns must be numeric, factor, or character coercible to numeric.")
  }

  out_list <- list()

  for (target_item in items) {
    other_items <- setdiff(items, target_item)
    rank_t <- to_rank_num(df[[target_item]])

    # ---------- Average rank ----------
    fit_avg <- fit_dim(rank_t)
    tb_avg  <- broom::tidy(fit_avg, conf.int = TRUE) |>
      dplyr::filter(.data$term == "D")

    b_avg <- -1 / (J - 1)
    avg_unit_est <- b_avg * tb_avg$estimate
    avg_unit_se  <- abs(b_avg) * tb_avg$std.error
    avg_unit_lo  <- pmin(b_avg * tb_avg$conf.low,  b_avg * tb_avg$conf.high)
    avg_unit_hi  <- pmax(b_avg * tb_avg$conf.low,  b_avg * tb_avg$conf.high)

    y_avg_unit <- (J - rank_t) / (J - 1)  # now definitely numeric in [0,1]
    mask_avg   <- (D == 0) & !is.na(y_avg_unit)

    sd_ctrl_avg_unit <- if (is.null(w)) {
      if (!any(mask_avg)) NA_real_ else {
        y0 <- y_avg_unit[mask_avg]
        n0 <- length(y0)
        if (n0 <= 1) NA_real_ else stats::sd(y0)
      }
    } else {
      wsd_masked_fsc(y_avg_unit, w_vec, mask_avg)
    }

    if (is.na(sd_ctrl_avg_unit) || sd_ctrl_avg_unit <= .Machine$double.eps) {
      est_std_avg <- se_std_avg <- lo_std_avg <- hi_std_avg <- NA_real_
    } else {
      est_std_avg <- avg_unit_est / sd_ctrl_avg_unit
      se_std_avg  <- avg_unit_se  / sd_ctrl_avg_unit
      lo_std_avg  <- avg_unit_lo  / sd_ctrl_avg_unit
      hi_std_avg  <- avg_unit_hi  / sd_ctrl_avg_unit
    }

    out_list[[length(out_list) + 1L]] <- tibble::tibble(
      estimate_raw   = tb_avg$estimate,
      std.error_raw  = tb_avg$std.error,
      conf.low_raw   = tb_avg$conf.low,
      conf.high_raw  = tb_avg$conf.high,
      estimate_unit  = avg_unit_est,
      std.error_unit = avg_unit_se,
      conf.low_unit  = avg_unit_lo,
      conf.high_unit = avg_unit_hi,
      estimate_std   = est_std_avg,
      std.error_std  = se_std_avg,
      conf.low_std   = lo_std_avg,
      conf.high_std  = hi_std_avg,
      sd_control_unit = sd_ctrl_avg_unit,
      p.value        = tb_avg$p.value,
      outcome        = paste0("Avg: ", target_item),
      effect_type    = "average rank",
      item           = target_item,
      J              = J
    )

    # ---------- Pairwise ----------
    for (oi in other_items) {
      y_pair <- as.integer(rank_t < df[[oi]])
      fit_pw <- fit_dim(y_pair)
      tb_pw  <- broom::tidy(fit_pw, conf.int = TRUE) |>
        dplyr::filter(.data$term == "D")

      est_unit <- tb_pw$estimate
      se_unit  <- tb_pw$std.error
      lo_unit  <- tb_pw$conf.low
      hi_unit  <- tb_pw$conf.high

      mask_pw <- ctrl_mask & !is.na(y_pair)
      p0 <- wmean_masked(y_pair, w_vec, mask_pw)
      sd0 <- if (is.na(p0)) NA_real_ else sqrt(pmax(p0 * (1 - p0), 0))

      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- est_unit / sd0
        se_std  <- se_unit  / sd0
        lo_std  <- lo_unit  / sd0
        hi_std  <- hi_unit  / sd0
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw   = tb_pw$estimate,
        std.error_raw  = tb_pw$std.error,
        conf.low_raw   = tb_pw$conf.low,
        conf.high_raw  = tb_pw$conf.high,
        estimate_unit  = est_unit,
        std.error_unit = se_unit,
        conf.low_unit  = lo_unit,
        conf.high_unit = hi_unit,
        estimate_std   = est_std,
        std.error_std  = se_std,
        conf.low_std   = lo_std,
        conf.high_std  = hi_std,
        sd_control_unit = sd0,
        p.value        = tb_pw$p.value,
        outcome        = paste0("v. ", oi),
        effect_type    = "pairwise ranking",
        item           = target_item,
        J              = J
      )
    }

    # ---------- Top-k ----------
    for (k in seq_len(top_k_max)) {
      y_topk <- as.integer(rank_t <= k)
      fit_tk <- fit_dim(y_topk)
      tb_tk  <- broom::tidy(fit_tk, conf.int = TRUE) |>
        dplyr::filter(.data$term == "D")

      est_unit <- tb_tk$estimate
      se_unit  <- tb_tk$std.error
      lo_unit  <- tb_tk$conf.low
      hi_unit  <- tb_tk$conf.high

      mask_tk <- ctrl_mask & !is.na(y_topk)
      p0 <- wmean_masked(y_topk, w_vec, mask_tk)
      sd0 <- if (is.na(p0)) NA_real_ else sqrt(pmax(p0 * (1 - p0), 0))

      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- est_unit / sd0
        se_std  <- se_unit  / sd0
        lo_std  <- lo_unit  / sd0
        hi_std  <- hi_unit  / sd0
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw   = tb_tk$estimate,
        std.error_raw  = tb_tk$std.error,
        conf.low_raw   = tb_tk$conf.low,
        conf.high_raw  = tb_tk$conf.high,
        estimate_unit  = est_unit,
        std.error_unit = se_unit,
        conf.low_unit  = lo_unit,
        conf.high_unit = hi_unit,
        estimate_std   = est_std,
        std.error_std  = se_std,
        conf.low_std   = lo_std,
        conf.high_std  = hi_std,
        sd_control_unit = sd0,
        p.value        = tb_tk$p.value,
        outcome        = paste0("Top-", k),
        effect_type    = "top-k ranking",
        item           = target_item,
        J              = J
      )
    }

    # ---------- Marginal ----------
    for (r in seq_len(J)) {
      y_marg <- as.integer(rank_t == r)
      fit_mg <- fit_dim(y_marg)
      tb_mg  <- broom::tidy(fit_mg, conf.int = TRUE) |>
        dplyr::filter(.data$term == "D")

      est_unit <- tb_mg$estimate
      se_unit  <- tb_mg$std.error
      lo_unit  <- tb_mg$conf.low
      hi_unit  <- tb_mg$conf.high

      mask_mg <- ctrl_mask & !is.na(y_marg)
      p0 <- wmean_masked(y_marg, w_vec, mask_mg)
      sd0 <- if (is.na(p0)) NA_real_ else sqrt(pmax(p0 * (1 - p0), 0))

      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- est_unit / sd0
        se_std  <- se_unit  / sd0
        lo_std  <- lo_unit  / sd0
        hi_std  <- hi_unit  / sd0
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw   = tb_mg$estimate,
        std.error_raw  = tb_mg$std.error,
        conf.low_raw   = tb_mg$conf.low,
        conf.high_raw  = tb_mg$conf.high,
        estimate_unit  = est_unit,
        std.error_unit = se_unit,
        conf.low_unit  = lo_unit,
        conf.high_unit = hi_unit,
        estimate_std   = est_std,
        std.error_std  = se_std,
        conf.low_std   = lo_std,
        conf.high_std  = hi_std,
        sd_control_unit = sd0,
        p.value        = tb_mg$p.value,
        outcome        = paste0("Ranked ", r),
        effect_type    = "marginal ranking",
        item           = target_item,
        J              = J
      )
    }
  }

  results_wide <- dplyr::bind_rows(out_list) |>
    dplyr::mutate(
      n             = n,
      n_treat       = n_treat,
      n_control     = n_ctrl,
      w_sum         = w_sum,
      w_treat_sum   = w_treat_sum,
      w_control_sum = w_ctrl_sum
    )

  results <- tidyr::pivot_longer(
    results_wide,
    cols = c(
      estimate_raw, estimate_unit, estimate_std,
      std.error_raw, std.error_unit, std.error_std,
      conf.low_raw, conf.low_unit, conf.low_std,
      conf.high_raw, conf.high_unit, conf.high_std
    ),
    names_to  = c(".value", "scale"),
    names_pattern = "(estimate|std\\.error|conf\\.low|conf\\.high)_(raw|unit|std)"
  ) |>
    dplyr::arrange(.data$item, .data$effect_type, .data$outcome, .data$scale) |>
    tibble::as_tibble()

  if (!return_models) return(results)
  list(results = results, models = NULL)
}
