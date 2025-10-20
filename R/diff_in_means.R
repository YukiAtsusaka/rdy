#' Estimate ranking effects (unweighted): raw, unit, and standardized (raw/sd0)
#'
#' For each item in `items`, constructs four outcome families:
#' 1) Average rank of the target (1 = best):
#'    - raw outcome: R in {1,...,J} (lower is better)
#'    - unit outcome: U = (J - R)/(J - 1) in [0,1] (higher is better)
#'    - std outcome: raw effect divided by SD(R | control)
#' 2) Pairwise wins vs. others (binary): Y = 1{R_target < R_other}
#' 3) Top-k (binary): Y = 1{R_target <= k}
#' 4) Marginal (binary): Y = 1{R_target == r}
#'
#' Estimation is unweighted difference-in-means via estimatr::lm_robust(y ~ D).
#' For each effect we report:
#'   - Raw ATE (for average-rank only; binaries have NA in raw block by design)
#'   - Unit ATE (for all; binaries' unit == raw)
#'   - Standardized ATE: (raw effect) / sd0_raw, where sd0_raw = SD of raw outcome among controls
#'
#' @param data A data.frame/tibble with `items` and `treat`.
#' @param items Character vector of ranking columns (1 = top/best). Length >= 2.
#' @param treat String; binary treatment column (0/1, logical, or factor codable to 0/1).
#' @param se_type SE type for estimatr::lm_robust (default "HC2").
#' @param clusters Optional string for clustering variable name in `data`.
#' @param na_action "drop" (default) or "none". If "drop", listwise on used cols.
#' @param top_k_max Optional integer (default J-1).
#' @param ... Extra args passed to estimatr::lm_robust.
#'
#' @return Tibble with one row per (item × effect). Columns:
#'   estimate_*, std.error_*, conf.low_*, conf.high_* for scales raw/unit/std,
#'   sd_control_raw, sd_control_unit, p.value, outcome, effect_type, item, J, n, n_treat, n_control.
#'
#' @importFrom rlang .data
#' @export
diff_in_means <- function(data,
                          items,
                          treat,
                          se_type = "HC2",
                          clusters = NULL,
                          na_action = c("drop", "none"),
                          top_k_max = NULL,
                          ...) {

  requireNamespace("estimatr", quietly = TRUE)
  requireNamespace("broom",    quietly = TRUE)
  requireNamespace("dplyr",    quietly = TRUE)
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

  used_cols <- unique(c(items, treat, clusters))
  df <- data[, used_cols, drop = FALSE]
  if (na_action == "drop") {
    keep <- stats::complete.cases(df)
    df <- df[keep, , drop = FALSE]
  }

  # Treatment coding to {0,1}
  D <- df[[treat]]
  if (is.logical(D)) D <- as.integer(D)
  if (is.factor(D))  D <- as.integer(D) - 1L
  if (!all(D %in% c(0, 1)))
    stop("`treat` must be binary/codable to {0,1} (0 = control, 1 = treated).")

  # Diagnostics
  n       <- nrow(df)
  n_treat <- sum(D == 1)
  n_ctrl  <- sum(D == 0)
  ctrl_mask <- (D == 0)

  # DIM helper
  fit_dim <- function(y) {
    estimatr::lm_robust(
      y ~ D,
      data     = df,
      se_type  = se_type,
      clusters = if (!is.null(clusters)) df[[clusters]] else NULL,
      ...
    )
  }

  # Safe coercion to numeric R in {1,...,J}
  to_rank_num <- function(x) {
    if (is.numeric(x))   return(as.numeric(x))
    if (is.factor(x))    return(suppressWarnings(as.numeric(as.character(x))))
    if (is.character(x)) return(suppressWarnings(as.numeric(x)))
    stop("Rank columns must be numeric, factor, or character coercible to numeric.")
  }

  # Unweighted SD among controls
  sd_control <- function(y, mask_ctrl) {
    y0 <- y[mask_ctrl & !is.na(y)]
    if (length(y0) <= 1) return(NA_real_)
    stats::sd(y0)
  }

  J <- length(items)
  if (is.null(top_k_max)) top_k_max <- J - 1L
  top_k_max <- max(1L, min(as.integer(top_k_max), J - 1L))

  out_list <- list()

  for (target_item in items) {
    other_items <- setdiff(items, target_item)

    # ----------------------------
    # 1) Average-rank outcomes
    # ----------------------------
    R <- to_rank_num(df[[target_item]])         # raw ranks {1,...,J}
    U <- (J - R) / (J - 1)                      # unit in [0,1]

    # raw effect (R)
    fit_avg_raw <- fit_dim(R)
    tb_avg_raw  <- broom::tidy(fit_avg_raw, conf.int = TRUE) |>
      dplyr::filter(.data$term == "D")

    # unit effect (U)
    fit_avg_unit <- fit_dim(U)
    tb_avg_unit  <- broom::tidy(fit_avg_unit, conf.int = TRUE) |>
      dplyr::filter(.data$term == "D")

    # control SDs (raw & unit)
    sd0_avg_raw  <- sd_control(R, ctrl_mask)
    sd0_avg_unit <- sd_control(U, ctrl_mask)

    # standardized: **raw effect divided by sd0_raw**
    if (is.na(sd0_avg_raw) || sd0_avg_raw <= .Machine$double.eps) {
      est_std_avg <- se_std_avg <- lo_std_avg <- hi_std_avg <- NA_real_
    } else {
      est_std_avg <- tb_avg_raw$estimate / sd0_avg_raw
      se_std_avg  <- tb_avg_raw$std.error / sd0_avg_raw
      lo_std_avg  <- tb_avg_raw$conf.low / sd0_avg_raw
      hi_std_avg  <- tb_avg_raw$conf.high / sd0_avg_raw
    }

    out_list[[length(out_list) + 1L]] <- tibble::tibble(
      estimate_raw    = tb_avg_raw$estimate,
      std.error_raw   = tb_avg_raw$std.error,
      conf.low_raw    = tb_avg_raw$conf.low,
      conf.high_raw   = tb_avg_raw$conf.high,
      estimate_unit   = tb_avg_unit$estimate,
      std.error_unit  = tb_avg_unit$std.error,
      conf.low_unit   = tb_avg_unit$conf.low,
      conf.high_unit  = tb_avg_unit$conf.high,
      estimate_std    = est_std_avg,
      std.error_std   = se_std_avg,
      conf.low_std    = lo_std_avg,
      conf.high_std   = hi_std_avg,
      sd_control_raw  = sd0_avg_raw,
      sd_control_unit = sd0_avg_unit,
      p.value         = tb_avg_unit$p.value, # (p-value for unit-scale effect; raw p-value can be added if needed)
      outcome         = paste0("Avg: ", target_item),
      effect_type     = "average rank",
      item            = target_item,
      J               = J,
      n               = n,
      n_treat         = n_treat,
      n_control       = n_ctrl
    )

    # ----------------------------
    # 2) Pairwise (binary) outcomes
    # ----------------------------
    for (oi in other_items) {
      R_oi <- to_rank_num(df[[oi]])
      Y_pw <- as.integer(R < R_oi)   # raw = unit for binary

      fit_pw <- fit_dim(Y_pw)
      tb_pw  <- broom::tidy(fit_pw, conf.int = TRUE) |>
        dplyr::filter(.data$term == "D")

      sd0_pw_raw  <- sd_control(Y_pw, ctrl_mask)
      sd0_pw_unit <- sd0_pw_raw

      if (is.na(sd0_pw_raw) || sd0_pw_raw <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- tb_pw$estimate / sd0_pw_raw
        se_std  <- tb_pw$std.error / sd0_pw_raw
        lo_std  <- tb_pw$conf.low / sd0_pw_raw
        hi_std  <- tb_pw$conf.high / sd0_pw_raw
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw    = NA_real_,  # no separate raw model for binary; raw==unit
        std.error_raw   = NA_real_,
        conf.low_raw    = NA_real_,
        conf.high_raw   = NA_real_,
        estimate_unit   = tb_pw$estimate,
        std.error_unit  = tb_pw$std.error,
        conf.low_unit   = tb_pw$conf.low,
        conf.high_unit  = tb_pw$conf.high,
        estimate_std    = est_std,
        std.error_std   = se_std,
        conf.low_std    = lo_std,
        conf.high_std   = hi_std,
        sd_control_raw  = sd0_pw_raw,
        sd_control_unit = sd0_pw_unit,
        p.value         = tb_pw$p.value,
        outcome         = paste0("v. ", oi),
        effect_type     = "pairwise ranking",
        item            = target_item,
        J               = J,
        n               = n,
        n_treat         = n_treat,
        n_control       = n_ctrl
      )
    }

    # ----------------------------
    # 3) Top-k (binary) outcomes
    # ----------------------------
    for (k in seq_len(top_k_max)) {
      Y_tk <- as.integer(R <= k)

      fit_tk <- fit_dim(Y_tk)
      tb_tk  <- broom::tidy(fit_tk, conf.int = TRUE) |>
        dplyr::filter(.data$term == "D")

      sd0_tk_raw  <- sd_control(Y_tk, ctrl_mask)
      sd0_tk_unit <- sd0_tk_raw

      if (is.na(sd0_tk_raw) || sd0_tk_raw <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- tb_tk$estimate / sd0_tk_raw
        se_std  <- tb_tk$std.error / sd0_tk_raw
        lo_std  <- tb_tk$conf.low / sd0_tk_raw
        hi_std  <- tb_tk$conf.high / sd0_tk_raw
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw    = NA_real_,
        std.error_raw   = NA_real_,
        conf.low_raw    = NA_real_,
        conf.high_raw   = NA_real_,
        estimate_unit   = tb_tk$estimate,
        std.error_unit  = tb_tk$std.error,
        conf.low_unit   = tb_tk$conf.low,
        conf.high_unit  = tb_tk$conf.high,
        estimate_std    = est_std,
        std.error_std   = se_std,
        conf.low_std    = lo_std,
        conf.high_std   = hi_std,
        sd_control_raw  = sd0_tk_raw,
        sd_control_unit = sd0_tk_unit,
        p.value         = tb_tk$p.value,
        outcome         = paste0("Top-", k),
        effect_type     = "top-k ranking",
        item            = target_item,
        J               = J,
        n               = n,
        n_treat         = n_treat,
        n_control       = n_ctrl
      )
    }

    # ----------------------------
    # 4) Marginal (binary) outcomes
    # ----------------------------
    for (r in seq_len(J)) {
      Y_mg <- as.integer(R == r)

      fit_mg <- fit_dim(Y_mg)
      tb_mg  <- broom::tidy(fit_mg, conf.int = TRUE) |>
        dplyr::filter(.data$term == "D")

      sd0_mg_raw  <- sd_control(Y_mg, ctrl_mask)
      sd0_mg_unit <- sd0_mg_raw

      if (is.na(sd0_mg_raw) || sd0_mg_raw <= .Machine$double.eps) {
        est_std <- se_std <- lo_std <- hi_std <- NA_real_
      } else {
        est_std <- tb_mg$estimate / sd0_mg_raw
        se_std  <- tb_mg$std.error / sd0_mg_raw
        lo_std  <- tb_mg$conf.low / sd0_mg_raw
        hi_std  <- tb_mg$conf.high / sd0_mg_raw
      }

      out_list[[length(out_list) + 1L]] <- tibble::tibble(
        estimate_raw    = NA_real_,
        std.error_raw   = NA_real_,
        conf.low_raw    = NA_real_,
        conf.high_raw   = NA_real_,
        estimate_unit   = tb_mg$estimate,
        std.error_unit  = tb_mg$std.error,
        conf.low_unit   = tb_mg$conf.low,
        conf.high_unit  = tb_mg$conf.high,
        estimate_std    = est_std,
        std.error_std   = se_std,
        conf.low_std    = lo_std,
        conf.high_std   = hi_std,
        sd_control_raw  = sd0_mg_raw,
        sd_control_unit = sd0_mg_unit,
        p.value         = tb_mg$p.value,
        outcome         = paste0("Ranked ", r),
        effect_type     = "marginal ranking",
        item            = target_item,
        J               = J,
        n               = n,
        n_treat         = n_treat,
        n_control       = n_ctrl
      )
    }
  }

  results_wide <- dplyr::bind_rows(out_list)

  # Long format (raw/unit/std blocks)
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
