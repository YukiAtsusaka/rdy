#' Visualize effects from \code{diff_in_means()}
#'
#' Produces a dot-and-interval plot of treatment effects, with configurable layout.
#'
#' @param estimates Tibble returned by \code{diff_in_means()} (long form) with
#'   columns including: \code{estimate}, \code{std.error}, \code{conf.low},
#'   \code{conf.high}, \code{scale}, \code{outcome}, \code{effect_type}, \code{item}, etc.
#' @param scale Which scale to plot: \code{"unit"} (default), \code{"raw"}, or \code{"std"}.
#' @param layout \code{"facets"} (default) to facet by effect type and item, or
#'   \code{"single_panel"} for one panel per item with separators between effect blocks.
#' @param palette Named colors for \code{c("Negative","Not_significant","Positive")}.
#' @param point_size Point size (default 2).
#' @param line_width Line width for CI bars (default 0.5).
#' @details
#' Average rank and pairwise ranking labels are alphabetically ordered (after
#' stripping prefixes), top-k and marginal are ordered numerically.
#' @importFrom rlang .data
#' @export
viz_diff <- function(estimates,
                      scale = c("unit", "raw", "std"),
                      layout = c("facets", "single_panel"),
                      palette = c(Negative = "#b0015a",
                                  Not_significant = "gray70",
                                  Positive = "#128ba0"),
                      point_size = 2,
                      line_width = 0.5) {

  requireNamespace("dplyr",   quietly = TRUE)
  requireNamespace("forcats", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("rlang",   quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)

  scale  <- match.arg(scale)
  layout <- match.arg(layout)

  # filter to one scale + significance label
  df <- estimates %>%
    dplyr::filter(.data$scale == !!scale) %>%
    dplyr::mutate(
      sig_label = dplyr::case_when(
        .data$conf.low  > 0 ~ "Positive",
        .data$conf.high < 0 ~ "Negative",
        TRUE               ~ "Not_significant"
      )
    )

  # canonical facet order
  effect_order <- c("average rank", "pairwise ranking", "top-k ranking", "marginal ranking")
  df <- df %>% dplyr::mutate(effect_type = factor(.data$effect_type, levels = effect_order))

  # --------- GLOBAL, PER-TYPE ORDERS (no per-item variation) ---------
  # 1) Average & Pairwise: alphabetical by cleaned label
  avgpair_order <- df %>%
    dplyr::filter(.data$effect_type %in% c("average rank","pairwise ranking")) %>%
    dplyr::distinct(effect_type, outcome) %>%
    dplyr::mutate(
      sort_key = dplyr::case_when(
        effect_type == "average rank"    ~ outcome %>%
          stringr::str_remove("^\\s*Avg:\\s*")  %>% stringr::str_squish() %>% stringr::str_to_lower(),
        effect_type == "pairwise ranking"~ outcome %>%
          stringr::str_remove("^\\s*v\\.?\\s*") %>% stringr::str_squish() %>% stringr::str_to_lower(),
        TRUE ~ outcome %>% stringr::str_squish() %>% stringr::str_to_lower()
      )
    ) %>%
    dplyr::group_by(effect_type) %>%
    dplyr::arrange(sort_key, .by_group = TRUE) %>%
    dplyr::mutate(order_in_block = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(effect_type, outcome, order_in_block)

  # 2) Top-k: numeric order
  topk_order <- df %>%
    dplyr::filter(.data$effect_type == "top-k ranking") %>%
    dplyr::distinct(effect_type, outcome) %>%
    dplyr::mutate(k_num = suppressWarnings(as.integer(stringr::str_extract(outcome, "\\d+")))) %>%
    dplyr::arrange(k_num) %>%
    dplyr::mutate(order_in_block = dplyr::row_number()) %>%
    dplyr::select(effect_type, outcome, order_in_block)

  # 3) Marginal: numeric order
  marg_order <- df %>%
    dplyr::filter(.data$effect_type == "marginal ranking") %>%
    dplyr::distinct(effect_type, outcome) %>%
    dplyr::mutate(rank_num = suppressWarnings(as.integer(stringr::str_extract(outcome, "\\d+")))) %>%
    dplyr::arrange(rank_num) %>%
    dplyr::mutate(order_in_block = dplyr::row_number()) %>%
    dplyr::select(effect_type, outcome, order_in_block)

  order_map <- dplyr::bind_rows(avgpair_order, topk_order, marg_order)

  # attach global per-type order and create an ordered factor
  df <- df %>%
    dplyr::left_join(order_map, by = c("effect_type","outcome")) %>%
    dplyr::mutate(outcome_ord = forcats::fct_reorder(.data$outcome, .data$order_in_block, .desc = FALSE))

  # ------------------- PLOT -------------------
  if (layout == "facets") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$outcome_ord, y = .data$estimate)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high, color = .data$sig_label),
        linewidth = line_width
      ) +
      ggplot2::geom_point(ggplot2::aes(color = .data$sig_label), size = point_size) +
      ggplot2::scale_color_manual(values = palette) +
      ggplot2::labs(x = NULL, y = NULL, color = NULL) +
      ggplot2::facet_grid(effect_type ~ item, scales = "free_y", space = "free_y") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        strip.text = ggplot2::element_text(size = 9),
        axis.text.y  = ggplot2::element_text(size = 8),
        axis.text.x  = ggplot2::element_text(size = 8)
      )
    return(p)
  }

  # single-panel with separators (unchanged; uses the same global order)
  df_sp <- df %>%
    dplyr::mutate(effect_rank = as.integer(.data$effect_type)) %>%
    dplyr::group_by(.data$item) %>%
    dplyr::arrange(.data$effect_rank, .data$order_in_block, .by_group = TRUE) %>%
    dplyr::mutate(outcome_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(outcome_blocked = forcats::fct_reorder(.data$outcome, .data$outcome_order, .desc = FALSE))

  sep_df <- df_sp %>%
    dplyr::distinct(item, outcome_blocked, effect_type, effect_rank) %>%
    dplyr::arrange(item, effect_rank, outcome_blocked) %>%
    dplyr::group_by(item, effect_rank) %>%
    dplyr::summarise(n_block = dplyr::n(), .groups = "drop_last") %>%
    dplyr::summarise(xints = cumsum(.data$n_block) + 0.5, .groups = "drop") %>%
    dplyr::group_by(item) %>%
    dplyr::mutate(xint = dplyr::lag(.data$xints)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$xint)) %>%
    dplyr::select(item, xint)

  p <- ggplot2::ggplot(df_sp, ggplot2::aes(x = .data$outcome_blocked, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high, color = .data$sig_label),
      linewidth = line_width
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data$sig_label), size = point_size) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::labs(x = NULL, y = NULL, color = NULL) +
    ggplot2::facet_grid(~ item, scales = "free_y", space = "free_y") +
    ggplot2::geom_vline(
      data = sep_df,
      ggplot2::aes(xintercept = .data$xint),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 9),
      axis.text.y  = ggplot2::element_text(size = 8),
      axis.text.x  = ggplot2::element_text(size = 8)
    )

  return(p)
}
