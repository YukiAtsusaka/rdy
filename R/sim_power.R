#' Monte Carlo Power Analysis for Difference-in-Means (with 95% CI)
#'
#' Repeatedly samples from a finite population, estimates difference-in-means
#' for outcomes matched by `items_pattern`, and computes power against the
#' population "truth" (taken as the point estimates from `diff_in_means(pop, ...)`).
#' Power is defined as: 95% CI excludes 0 AND the estimate's sign matches the truth.
#'
#' @param pop A data.frame (finite population).
#' @param items_pattern Regex to select outcome columns in `pop`. Default "^ch_".
#' @param treat Name of the treatment indicator column in `pop`. Default "treat".
#' @param n_sample Sample size per replication. Default 500.
#' @param n_rep Number of Monte Carlo replications. Default 50.
#' @param scale_filter Character scalar; if not NULL, filter truth/estimates to this scale (e.g., "std").
#' @param eps Numerical tolerance for defining "true_nonzero". Default 0.
#' @param seed Optional integer seed for reproducibility. Default NULL.
#' @return A list with:
#' \itemize{
#'   \item \code{sim_results}: tibble of per-replication indicators and metadata
#'   \item \code{summary}: tibble of power by (item, outcome, effect_type)
#'   \item \code{plot}: ggplot object showing power by outcome, faceted by item
#' }
#' @details
#' This function expects that \code{diff_in_means()} returns a tibble including
#' columns: \code{item}, \code{outcome}, \code{effect_type}, \code{scale}, \code{estimate},
#' \code{conf.low}, \code{conf.high}. Missing estimates/intervals are handled gracefully.
#'
#' @examples
#' \dontrun{
#' res <- sim_power(pop, items_pattern = "^ch_", treat = "treat",
#'                      n_sample = 500, n_rep = 50, scale_filter = "std", seed = 123)
#' res$summary
#' print(res$plot)
#' }
#' @export
sim_power <- function(pop,
                      items_pattern = "^ch_",
                      treat = "treat",
                      n_sample = 500,
                      n_rep = 50,
                      scale_filter = "std",
                      eps = 0,
                      seed = NULL) {

  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("purrr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("readr", quietly = TRUE)   # if you load pop elsewhere, this isn't needed
  requireNamespace("scales", quietly = TRUE)

  if (!is.null(seed)) set.seed(seed)

  # Identify item columns
  items <- grep(items_pattern, colnames(pop), value = TRUE)
  if (length(items) == 0L) {
    stop("No outcome columns matched by `items_pattern` in `pop`.")
  }
  if (!treat %in% names(pop)) {
    stop("`treat` column not found in `pop`.")
  }

  # Population "truth"
  out <- diff_in_means(pop, items, treat = treat)
  if (!is.null(scale_filter)) {
    if (!"scale" %in% names(out)) {
      stop("`scale_filter` provided but `diff_in_means()` output has no `scale` column.")
    }
    out <- dplyr::filter(out, .data$scale == scale_filter)
  }

  # Keep only necessary truth columns
  truth <- out %>%
    dplyr::select(dplyr::any_of(c("item", "outcome", "effect_type", "scale", "estimate"))) %>%
    dplyr::rename(true = "estimate")

  # Simulation loop
  sim_results <- purrr::map_dfr(seq_len(n_rep), function(r) {
    samp <- dplyr::slice_sample(pop, n = n_sample)

    est <- diff_in_means(samp, items, treat = treat)
    if (!is.null(scale_filter)) {
      est <- dplyr::filter(est, .data$scale == scale_filter)
    }

    cmp <- est %>%
      dplyr::select(dplyr::any_of(c("item", "outcome", "effect_type", "scale",
                                    "estimate", "conf.low", "conf.high"))) %>%
      dplyr::left_join(truth, by = c("item", "outcome", "effect_type", "scale")) %>%
      dplyr::mutate(
        true_nonzero = dplyr::if_else(!is.na(.data$true), abs(.data$true) > eps, NA),
        sig = dplyr::if_else(
          !is.na(.data$conf.low) & !is.na(.data$conf.high),
          (.data$conf.low > 0 | .data$conf.high < 0),
          NA
        ),
        correct_sign = dplyr::if_else(
          !is.na(.data$estimate) & !is.na(.data$true),
          sign(.data$estimate) == sign(.data$true),
          NA
        ),
        power_success = dplyr::if_else(
          !is.na(.data$sig) & !is.na(.data$correct_sign),
          .data$sig & .data$correct_sign,
          NA
        ),
        rep = r
      )

    cmp
  })

  # Aggregate across simulations
  finite_summary <- sim_results %>%
    dplyr::group_by(.data$item, .data$outcome, .data$effect_type) %>%
    dplyr::summarise(
      power = mean(.data$power_success, na.rm = TRUE),
      .groups = "drop"
    )

  # Plot
  plt <- finite_summary %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$outcome, y = .data$power, fill = .data$item, alpha = .data$power)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::facet_wrap(~ .data$item) +
    ggplot2::labs(x = NULL, y = "Statistical Power") +
    ggplot2::scale_fill_manual(values = c("darkcyan", "maroon", "gray50", "gold4")) +
    ggplot2::scale_alpha(range = c(0.3, 1)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::ylim(0, 1) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.3,
                        color = scales::alpha("black", 0.5)) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::coord_flip()

  list(sim_results = sim_results, summary = finite_summary, plot = plt)
}
