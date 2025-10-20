#' Estimate ranking effects (UNWEIGHTED): precompute outcomes → fit → summarize
#' Returns: list(results = <tibble>, outcomes_long = <tibble>)
diff_in_means <- function(data,
                          items,
                          treat,
                          se_type = "HC2",
                          clusters = NULL,
                          na_action = c("none", "drop"),
                          top_k_max = NULL,
                          ...) {
  requireNamespace("dplyr",    quietly = TRUE)
  requireNamespace("tibble",   quietly = TRUE)
  requireNamespace("tidyr",    quietly = TRUE)
  requireNamespace("purrr",    quietly = TRUE)
  requireNamespace("broom",    quietly = TRUE)
  requireNamespace("estimatr", quietly = TRUE)

  na_action <- match.arg(na_action)
  stopifnot(is.data.frame(data))

  # ---- Checks ----
  if (!is.character(items) || length(items) < 2L)
    stop("`items` must be a character vector with at least two ranking columns.")
  if (!is.character(treat) || length(treat) != 1L)
    stop("`treat` must be a single column name (string).")
  if (!all(items %in% names(data))) stop("Some `items` not found in `data`.")
  if (!treat %in% names(data))      stop("`treat` not found in `data`.")
  if (!is.null(clusters) && !clusters %in% names(data))
    stop("`clusters` not found in `data`.")

  needed <- unique(c(items, treat, clusters))
  df <- data[, needed, drop = FALSE]

  # Up-front drop only on treat/clusters; outcomes handle their own NAs
  if (na_action == "drop") {
    keep <- stats::complete.cases(df[, c(treat, clusters), drop = FALSE])
    df <- df[keep, , drop = FALSE]
  }

  # Treat → {0,1}
  D <- df[[treat]]
  if (is.logical(D)) D <- as.integer(D)
  if (is.factor(D))  D <- as.integer(D) - 1L
  if (!all(D %in% c(0, 1), na.rm = TRUE))
    stop("`treat` must be binary/codable to {0,1} (0=control, 1=treated).")

  # Diagnostics (overall; per-outcome N may differ)
  n_all       <- nrow(df)
  n_treat_all <- sum(D == 1, na.rm = TRUE)
  n_ctrl_all  <- sum(D == 0, na.rm = TRUE)

  # Helpers --------------------------------------------------------------------
  to_rank_num <- function(x) {
    if (is.numeric(x))   return(as.numeric(x))
    if (is.factor(x))    return(suppressWarnings(as.numeric(as.character(x))))
    if (is.character(x)) return(suppressWarnings(as.numeric(x)))
    stop("Rank columns must be numeric, factor, or character coercible to numeric.")
  }

  # Always have a cluster vector; NA if no clustering
  Cl_vec <- if (!is.null(clusters)) df[[clusters]] else rep(NA, nrow(df))

  # Fit y ~ D and compute SDs using exactly the rows used in that regression
  fit_with_sds <- function(y_fit, y_raw_for_sd, Dvec, Cl = NULL) {
    use_cluster <- !is.null(Cl) && any(!is.na(Cl))
    rows_used <- (!is.na(y_fit)) & (!is.na(Dvec))
    if (use_cluster) rows_used <- rows_used & (!is.na(Cl))

    if (sum(rows_used) < 2L) {
      return(tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_,
                            p.value = NA_real_,
                            sd_control_raw = NA_real_,
                            sd_control_current = NA_real_))
    }

    y_ <- y_fit[rows_used]
    D_ <- Dvec[rows_used]
    ok_groups <- any(D_ == 0) && any(D_ == 1)

    # SDs among controls on the same rows
    sd0_raw     <- stats::sd(y_raw_for_sd[rows_used][D_ == 0], na.rm = TRUE)
    sd0_current <- stats::sd(y_[D_ == 0],                      na.rm = TRUE)

    if (!ok_groups) {
      return(tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_,
                            p.value = NA_real_,
                            sd_control_raw = sd0_raw,
                            sd_control_current = sd0_current))
    }

    if (!use_cluster) {
      dat <- data.frame(y = y_, D = D_)
      fit <- tryCatch(estimatr::lm_robust(y ~ D, data = dat, se_type = se_type, ...),
                      error = function(e) NULL)
    } else {
      Cl_ <- Cl[rows_used]
      dat <- data.frame(y = y_, D = D_, Cl = Cl_)
      fit <- tryCatch(estimatr::lm_robust(y ~ D, data = dat,
                                          se_type = se_type, clusters = dat$Cl, ...),
                      error = function(e) NULL)
    }

    if (is.null(fit)) {
      return(tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_,
                            p.value = NA_real_,
                            sd_control_raw = sd0_raw,
                            sd_control_current = sd0_current))
    }

    tb <- broom::tidy(fit, conf.int = TRUE)
    tb <- tb[tb$term == "D", c("estimate","std.error","conf.low","conf.high","p.value")]
    if (nrow(tb) == 0) {
      tb <- tibble::tibble(estimate = NA_real_, std.error = NA_real_,
                           conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
    }

    dplyr::bind_cols(tb, tibble::tibble(sd_control_raw = sd0_raw,
                                        sd_control_current = sd0_current))
  }

  J <- length(items)
  if (is.null(top_k_max)) top_k_max <- J - 1L
  top_k_max <- max(1L, min(as.integer(top_k_max), J - 1L))

  # ============= (1) PRECOMPUTE OUTCOMES (LONG) ===============================
  # For every effect, create BOTH raw and unit rows.
  # - Average rank: raw = R; unit = U = (J - R)/(J - 1)
  # - Binary families: raw == unit = Y (we duplicate rows for raw & unit)
  # y_raw_for_sd is the raw outcome we standardize by (R for avg-rank; Y for binaries)
  outcomes <- list()

  for (it in items) {
    R <- to_rank_num(df[[it]])
    U <- (J - R) / (J - 1)

    # Average rank
    outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
      item = it, effect_type = "average rank", outcome = paste0("Avg: ", it),
      scale = "raw",  y = R, y_raw_for_sd = R, D = D, Cl = Cl_vec, J = J
    )
    outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
      item = it, effect_type = "average rank", outcome = paste0("Avg: ", it),
      scale = "unit", y = U, y_raw_for_sd = R, D = D, Cl = Cl_vec, J = J
    )

    # Pairwise (binary): create raw and unit (identical)
    for (oi in setdiff(items, it)) {
      R_oi <- to_rank_num(df[[oi]])
      Y_pw <- as.integer(R < R_oi)
      # raw
      outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
        item = it, effect_type = "pairwise ranking", outcome = paste0("v. ", oi),
        scale = "raw", y = Y_pw, y_raw_for_sd = Y_pw, D = D, Cl = Cl_vec, J = J
      )
      # unit (duplicate)
      outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
        item = it, effect_type = "pairwise ranking", outcome = paste0("v. ", oi),
        scale = "unit", y = Y_pw, y_raw_for_sd = Y_pw, D = D, Cl = Cl_vec, J = J
      )
    }

    # Top-k (binary): raw + unit (identical)
    for (k in seq_len(top_k_max)) {
      Y_tk <- as.integer(R <= k)
      outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
        item = it, effect_type = "top-k ranking", outcome = paste0("Top-", k),
        scale = "raw", y = Y_tk, y_raw_for_sd = Y_tk, D = D, Cl = Cl_vec, J = J
      )
      outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
        item = it, effect_type = "top-k ranking", outcome = paste0("Top-", k),
        scale = "unit", y = Y_tk, y_raw_for_sd = Y_tk, D = D, Cl = Cl_vec, J = J
      )
    }

    # Marginal (binary): raw + unit (identical)
    for (r in seq_len(J)) {
      Y_mg <- as.integer(R == r)
      outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
        item = it, effect_type = "marginal ranking", outcome = paste0("Ranked ", r),
        scale = "raw", y = Y_mg, y_raw_for_sd = Y_mg, D = D, Cl = Cl_vec, J = J
      )
      outcomes[[length(outcomes) + 1L]] <- tibble::tibble(
        item = it, effect_type = "marginal ranking", outcome = paste0("Ranked ", r),
        scale = "unit", y = Y_mg, y_raw_for_sd = Y_mg, D = D, Cl = Cl_vec, J = J
      )
    }
  }

  outcomes_long <- dplyr::bind_rows(outcomes)

  # ============= (2) FIT PER OUTCOME (and compute SDs on same rows) ===========
  fitted <- outcomes_long |>
    dplyr::group_by(item, effect_type, outcome, scale) |>
    dplyr::summarise(
      stats = list(fit_with_sds(y_fit = y, y_raw_for_sd = y_raw_for_sd,
                                Dvec = D, Cl = Cl)),
      .groups = "drop"
    ) |>
    tidyr::unnest_wider(col = stats, names_sep = "_")

  # ============= (3) ADD STANDARDIZED ROWS ====================================
  fitted_core <- fitted |>
    dplyr::mutate(
      J = unique(outcomes_long$J),
      n_all = n_all, n_treat_all = n_treat_all, n_control_all = n_ctrl_all
    ) |>
    dplyr::rename(
      estimate        = stats_estimate,
      std.error       = stats_std.error,
      conf.low        = stats_conf.low,
      conf.high       = stats_conf.high,
      p.value         = stats_p.value,
      sd_control_raw  = stats_sd_control_raw,
      sd_control_unit = stats_sd_control_current  # for unit rows this equals the unit SD
    )

  # Standardize ALWAYS from the raw row of each effect (for binaries: raw==unit)
  base_raw <- fitted_core |> dplyr::filter(scale == "raw")

  mk_std <- function(base_df) {
    if (nrow(base_df) == 0) return(base_df[0, ])
    base_df |>
      dplyr::mutate(
        scale     = "std",
        estimate  = ifelse(is.na(sd_control_raw) | sd_control_raw <= .Machine$double.eps,
                           NA_real_, estimate / sd_control_raw),
        std.error = ifelse(is.na(sd_control_raw) | sd_control_raw <= .Machine$double.eps,
                           NA_real_, std.error / sd_control_raw),
        conf.low  = ifelse(is.na(sd_control_raw) | sd_control_raw <= .Machine$double.eps,
                           NA_real_, conf.low / sd_control_raw),
        conf.high = ifelse(is.na(sd_control_raw) | sd_control_raw <= .Machine$double.eps,
                           NA_real_, conf.high / sd_control_raw)
      )
  }

  std_rows <- mk_std(base_raw)

  results <- dplyr::bind_rows(
    fitted_core,
    std_rows
  ) |>
    dplyr::arrange(item, effect_type, outcome,
                   factor(scale, levels = c("raw","unit","std")))

  list(
    results = results,
    outcomes_long = outcomes_long
  )
}
