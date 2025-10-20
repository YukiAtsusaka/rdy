#' Estimate ranking effects (unweighted): raw, unit, and standardized (raw/sd0)
#'
#' For each item in `items`, constructs four outcome families:
#' 1) Average rank (1 = best):
#'    - raw: R in {1,...,J} (lower is better)
#'    - unit: U = (J - R)/(J - 1) in [0,1] (higher is better)
#'    - std: raw effect divided by SD(R | control)
#' 2) Pairwise wins vs. others (binary): Y = 1{R_target < R_other}
#' 3) Top-k (binary): Y = 1{R_target <= k}
#' 4) Marginal (binary): Y = 1{R_target == r}
#'
#' Estimation is unweighted DIM via estimatr::lm_robust(y ~ D), fit on
#' per-outcome subsets (no global listwise drop).
#'
#' @param data A data.frame/tibble with `items` and `treat`.
#' @param items Character vector of ranking columns (1 = top/best). Length >= 2.
#' @param treat String; binary treatment column (0/1, logical, or factor codable to 0/1).
#' @param se_type SE type for estimatr::lm_robust (default "HC2").
#' @param clusters Optional string for clustering variable name in `data`.
#' @param na_action "none" (default) or "drop". If "drop", only `treat`/`clusters` are filtered globally;
#'                  outcomes still do per-outcome NA filtering.
#' @param top_k_max Optional integer (default J-1).
#' @param ... Extra args passed to estimatr::lm_robust.
#'
#' @return Tibble with rows for all effects across all items and scales.
#'         Columns: estimate_*, std.error_*, conf.low_*, conf.high_* (raw/unit/std),
#'         sd_control_raw, sd_control_unit, p.value, outcome, effect_type, item, J,
#'         n, n_treat, n_control.
#'
#' @importFrom rlang .data
#' @export
diff_in_means <- function(data,
                          items,
                          treat,
                          se_type = "HC2",
                          clusters = NULL,
                          na_action = c("none", "drop"),
                          top_k_max = NULL,
                          ...) {

  requireNamespace("estimatr", quietly = TRUE)
  requireNamespace("broom",    quietly = TRUE)
  requireNamespace("dplyr",    quietly = TRUE)
  requireNamespace("tibble",   quietly = TRUE)
  requireNamespace("tidyr",    quietly = TRUE)

  na_action <- match.arg(na_action)
  stopifnot(is.data.frame(data))

  # ---- Basic checks -----------------------------------------------------------
  if (!is.character(items) || length(items) < 2L)
    stop("`items` must be a character vector with at least two ranking columns.")
  if (!is.character(treat) || length(treat) != 1L)
    stop("`treat` must be a single column name (string).")
  if (!all(items %in% names(data))) stop("Some `items` not found in `data`.")
  if (!treat %in% names(data))      stop("`treat` not found in `data`.")
  if (!is.null(clusters) && !clusters %in% names(data))
    stop("`clusters` not found in `data`.")

  used_cols <- unique(c(items, treat, clusters))
  df <- data[, used_cols, drop = FALSE]

  # Global filtering only for treat/clusters if requested; outcomes handle NAs themselves
  if (na_action == "drop") {
    keep <- stats::complete.cases(df[, c(treat, clusters), drop = FALSE])
    df <- df[keep, , drop = FALSE]
  }

  # Treatment to {0,1}
  D <- df[[treat]]
  if (is.logical(D)) D <- as.integer(D)
  if (is.factor(D))  D <- as.integer(D) - 1L
  if (!all(D %in% c(0, 1), na.rm = TRUE))
    stop("`treat` must be binary/codable to {0,1} (0 = control, 1 = treated).")

  # Overall diagnostics (note: per-outcome N may differ)
  n_all       <- nrow(df)
  n_treat_all <- sum(D == 1, na.rm = TRUE)
  n_ctrl_all  <- sum(D == 0, na.rm = TRUE)

  # Helpers --------------------------------------------------------------------
  # Safe coercion to numeric R in {1,...,J}
  to_rank_num <- function(x) {
    if (is.numeric(x))   return(as.numeric(x))
    if (is.factor(x))    return(suppressWarnings(as.numeric(as.character(x))))
    if (is.character(x)) return(suppressWarnings(as.numeric(x)))
    stop("Rank columns must be numeric, factor, or character coercible to numeric.")
  }

  # Unweighted SD among controls on *given* scale
  sd_control <- function(y, Dvec) {
    m <- (!is.na(y)) & (!is.na(Dvec)) & (Dvec == 0)
    y0 <- y[m]
    if (length(y0) <= 1) return(NA_real_)
    stats::sd(y0)
  }

  # Safe per-outcome DIM fit y ~ D with optional cluster
  safe_fit <- function(y, Dvec, clusters_vec = NULL) {
    # Per-outcome mask
    mask <- (!is.na(y)) & (!is.na(Dvec))
    if (!is.null(clusters_vec)) mask <- mask & (!is.na(clusters_vec))

    if (sum(mask) < 2L) {
      return(tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_,
                            p.value = NA_real_))
    }
    y_sub <- y[mask]
    D_sub <- Dvec[mask]
    # Need both groups present
    if (!(any(D_sub == 0) && any(D_sub == 1))) {
      return(tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_,
                            p.value = NA_real_))
    }

    if (is.null(clusters_vec)) {
      dat <- data.frame(y = y_sub, D = D_sub)
      fit <- tryCatch(
        estimatr::lm_robust(y ~ D, data = dat, se_type = se_type, ...),
        error = function(e) NULL
      )
    } else {
      Cl_sub <- clusters_vec[mask]
      dat <- data.frame(y = y_sub, D = D_sub, Cl = Cl_sub)
      fit <- tryCatch(
        estimatr::lm_robust(y ~ D, data = dat, se_type = se_type, clusters = Cl, ...),
        error = function(e) NULL
      )
    }

    if (is.null(fit)) {
      return(tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_,
                            p.value = NA_real_))
    }
    tb <- broom::tidy(fit, conf.int = TRUE)
    tb <- tb[tb$term == "D", c("estimate","std.error","conf.low","conf.high","p.value")]
    if (nrow(tb) == 0) {
      tb <- tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                           conf.low = NA_real_, conf.high = NA_real_,
                           p.value = NA_real_)
    }
    tibble::as_tibble(tb)
  }

  J <- length(items)
  if (is.null(top_k_max)) top_k_max <- J - 1L
  top_k_max <- max(1L, min(as.integer(top_k_max), J - 1L))

  Cl_vec <- if (!is.null(clusters)) df[[clusters]] else NULL

  out_list <- list()

  for (target_item in items) {
    other_items <- setdiff(items, target_item)

    # ----- Average rank: raw (R), unit (U), std = raw/sd0_raw ------------------
    R <- to_rank_num(df[[target_item]])
    U <- (J - R) / (J - 1)

    tb_raw  <- safe_fit(R, D, Cl_vec)
    tb_unit <- safe_fit(U, D, Cl_vec)

    sd0_raw  <- sd_control(R, D)
    sd0_unit <- sd_control(U, D)

    if (is.na(sd0_raw) || sd0_raw <= .Machine$double.eps) {
      est_std <- se_std <- lo_std <- hi_std <- NA_real_
    } else {
      est_std <- tb_raw$estimate / sd0_raw
      se_std  <- tb_raw$std.error / sd0_raw
      lo_std  <- tb_raw$conf.low / sd0_raw
      hi_std  <- tb_raw$conf.high / sd0_raw
    }

    out_list[[length(out_list) + 1L]] <- tibble::tibble(
      estimate_raw     = tb_raw$estimate,
      std.error_raw    = tb_raw$std.error,
      conf.low_raw     = tb_raw$conf.low,
      conf.high_raw    = tb_raw$conf.high,
      estimate_unit    = tb_unit$estimate,
      std.error_unit   = tb_unit$std.error,
      conf.low_unit    = tb_unit$conf.low,
      conf.high_unit   = tb_unit$conf.high,
      estimate_std     = est_std,
      std.error_std    = se_std,
      conf.low_std     = lo_std,
      conf.high_std    = hi_std,
      sd_control_raw   = sd0_raw,
      sd_control_unit  = sd0_unit,
      p.value          = tb_unit$p.value,  # p-value for unit-scale effect
      outcome          = paste0("Avg: ", target_item),
      effect_type      = "average rank",
      item             = target_item,
      J                = J,
      n                = n_all,
      n_treat          = n_treat_all,
      n_control        = n_ctrl_all
    )

    # ----- Pairwise (binary): raw==unit; std = raw/sd0_raw ---------------------
    for (oi in other_items) {
      R_oi <- to_rank_num(df[[oi]])
      Y_pw <- as.integer(R < R_oi)

      tb <- safe_fit(Y_pw, D, Cl_vec)
      sd0 <- sd_control(Y_pw, D)

      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- tb$estimate / sd0
        se_std  <- tb$std.error / sd0
        lo_std  <- tb$conf.low / sd0
        hi_std  <- tb$conf.high / sd0
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw     = NA_real_,  # raw model not separate for binary
        std.error_raw    = NA_real_,
        conf.low_raw     = NA_real_,
        conf.high_raw    = NA_real_,
        estimate_unit    = tb$estimate,
        std.error_unit   = tb$std.error,
        conf.low_unit    = tb$conf.low,
        conf.high_unit   = tb$conf.high,
        estimate_std     = est_std,
        std.error_std    = se_std,
        conf.low_std     = lo_std,
        conf.high_std    = hi_std,
        sd_control_raw   = sd0,
        sd_control_unit  = sd0,
        p.value          = tb$p.value,
        outcome          = paste0("v. ", oi),
        effect_type      = "pairwise ranking",
        item             = target_item,
        J                = J,
        n                = n_all,
        n_treat          = n_treat_all,
        n_control        = n_ctrl_all
      )
    }

    # ----- Top-k (binary) ------------------------------------------------------
    for (k in seq_len(top_k_max)) {
      Y_tk <- as.integer(R <= k)

      tb <- safe_fit(Y_tk, D, Cl_vec)
      sd0 <- sd_control(Y_tk, D)

      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- tb$estimate / sd0
        se_std  <- tb$std.error / sd0
        lo_std  <- tb$conf.low / sd0
        hi_std  <- tb$conf.high / sd0
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw     = NA_real_,
        std.error_raw    = NA_real_,
        conf.low_raw     = NA_real_,
        conf.high_raw    = NA_real_,
        estimate_unit    = tb$estimate,
        std.error_unit   = tb$std.error,
        conf.low_unit    = tb$conf.low,
        conf.high_unit   = tb$conf.high,
        estimate_std     = est_std,
        std.error_std    = se_std,
        conf.low_std     = lo_std,
        conf.high_std    = hi_std,
        sd_control_raw   = sd0,
        sd_control_unit  = sd0,
        p.value          = tb$p.value,
        outcome          = paste0("Top-", k),
        effect_type      = "top-k ranking",
        item             = target_item,
        J                = J,
        n                = n_all,
        n_treat          = n_treat_all,
        n_control        = n_ctrl_all
      )
    }

    # ----- Marginal (binary) ---------------------------------------------------
    for (r in seq_len(J)) {
      Y_mg <- as.integer(R == r)

      tb <- safe_fit(Y_mg, D, Cl_vec)
      sd0 <- sd_control(Y_mg, D)

      if (is.na(sd0) || sd0 <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- tb$estimate / sd0
        se_std  <- tb$std.error / sd0
        lo_std  <- tb$conf.low / sd0
        hi_std  <- tb$conf.high / sd0
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw     = NA_real_,
        std.error_raw    = NA_real_,
        conf.low_raw     = NA_real_,
        conf.high_raw    = NA_real_,
        estimate_unit    = tb$estimate,
        std.error_unit   = tb$std.error,
        conf.low_unit    = tb$conf.low,
        conf.high_unit   = tb$conf.high,
        estimate_std     = est_std,
        std.error_std    = se_std,
        conf.low_std     = lo_std,
        conf.high_std    = hi_std,
        sd_control_raw   = sd0,
        sd_control_unit  = sd0,
        p.value          = tb$p.value,
        outcome          = paste0("Ranked ", r),
        effect_type      = "marginal ranking",
        item             = target_item,
        J                = J,
        n                = n_all,
        n_treat          = n_treat_all,
        n_control        = n_ctrl_all
      )
    }
  }

  results_wide <- dplyr::bind_rows(out_list)

  # Long format (raw/unit/std)
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

  results
}
