#' Ordinal Helpers & Plackett–Luce Sampler
#'
#' @name rpluce
#' @keywords sampling ranking
NULL

# ------------------------------------------------------------------------
# ordinal_seq -------------------------------------------------------------
# ------------------------------------------------------------------------

#' Generate an Ordinal Sequence from a Number
#'
#' Returns \code{c("1st","2nd","3rd", ...)} of length \code{n}.
#'
#' @param n A single numeric value: the desired length of the sequence.
#' @return A character vector of ordinal strings of length \code{n}.
#'
#' @examples
#' ordinal_seq(11)
#'
#' @export
ordinal_seq <- function(n) {
  # Sanity check on `n` argument
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || n < 0) {
    stop("`n` must be a non-negative numeric value of length 1.")
  }
  n <- as.integer(n)

  ordinal_suffix <- function(k) {
    if (k %% 100L %in% c(11L, 12L, 13L)) {
      "th"
    } else {
      c("th", "st", "nd", "rd", rep("th", 6L))[(k %% 10L) + 1L]
    }
  }

  seq_vec <- seq.int(1L, length.out = n)
  paste0(seq_vec, vapply(seq_vec, ordinal_suffix, FUN.VALUE = character(1)))
}

# ------------------------------------------------------------------------
# unique_alphabets --------------------------------------------------------
# ------------------------------------------------------------------------

#' Internal: Create a Vector of Unique Item Labels
#'
#' Generates up to the requested length of labels using \code{letters}, then
#' cartesian concatenations (a, b, ..., z, aa, ab, ..., etc.).
#'
#' @param m Integer: number of labels to return.
#' @return Character vector of length \code{m}.
#' @noRd
unique_alphabets <- function(m) {
  if (length(m) != 1L || !is.numeric(m) || is.na(m) || m < 1L) {
    stop("`m` must be a positive numeric value of length 1.")
  }
  m <- as.integer(m)
  if (m > 100000L) {
    stop("Reconsider: do you really need more than 100,000 items?")
  }

  if (m <= 26L) {
    return(letters[seq_len(m)])
  }

  # Build labels by expanding letter prefixes until we have >= m
  out <- letters
  block <- letters
  while (length(out) < m) {
    block <- as.vector(outer(block, letters, paste0))
    out <- c(out, block)
  }
  out[seq_len(m)]
}

# ------------------------------------------------------------------------
# rpluce ------------------------------------------------------------------
# ------------------------------------------------------------------------

#' Draw Samples from the Plackett–Luce Model
#'
#' Implements Algorithm 2.1 ("Efficient Sampling from Plackett–Luce")
#' in Xia (2019), §2.2.3. Draws full rankings for \code{t} items with
#' choice probabilities \code{prob}.
#'
#' @param n Integer: number of samples (assessors).
#' @param t Integer: number of items/alternatives to rank (\eqn{t \ge 1}).
#' @param prob Numeric vector of choice probabilities of length \code{t}
#'   (nonnegative, sums to 1 within tolerance).
#' @param choices Optional character vector of item labels of length \code{t}.
#'   If \code{NULL}, uses \code{unique_alphabets(t)}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @importFrom stats rmultinom
#' @return A data frame with \code{n} rows and \code{t} columns. Each row is a
#'   complete ranking (1st through t-th). Column names are \code{ordinal_seq(t)}.
#'
#' @examples
#' rpluce(n = 10, t = 3, prob = c(0.5, 0.3, 0.2), seed = 123)
#'
#' @export
rpluce <- function(n, t, prob, choices = NULL, seed = NULL) {
  # Seed
  if (!is.null(seed)) set.seed(seed)

  # Basic checks
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || n < 1L) {
    stop("`n` must be a positive numeric value of length 1.")
  }
  n <- as.integer(n)

  if (length(t) != 1L || !is.numeric(t) || is.na(t) || t < 1L) {
    stop("`t` must be a positive numeric value of length 1.")
  }
  t <- as.integer(t)

  if (!is.numeric(prob) || any(is.na(prob)) || length(prob) != t || any(prob < 0)) {
    stop("`prob` must be a nonnegative numeric vector of length `t` with no NAs.")
  }
  if (!isTRUE(all.equal(sum(prob), 1, tolerance = 1e-8))) {
    stop("`prob` must sum to 1 (within a small numerical tolerance).")
  }

  if (!is.null(choices)) {
    if (!is.character(choices) || length(choices) != t) {
      stop("`choices` must be a character vector of length `t`.")
    }
    item_labels <- choices
  } else {
    item_labels <- unique_alphabets(t)
  }

  # Handle t == 1 as a degenerate single-item case
  if (t == 1L) {
    # Only one possible ranking; repeat it n times
    res <- data.frame(`1st` = rep.int(item_labels[1L], n), check.names = FALSE)
    return(res)
  }

  # Storage
  rank_matrix <- matrix(NA_character_, nrow = n, ncol = t)

  # Sampling per assessor
  for (j in seq_len(n)) {
    R <- character(t)   # will fill 1..t

    # Initial choice set & probabilities
    A <- item_labels
    Gamma <- prob

    # Fill positions 1..(t-1)
    for (i in seq_len(t - 1L)) {
      draw <- as.vector(stats::rmultinom(n = 1L, size = 1L, prob = Gamma))
      chosen_idx <- which.max(draw)     # index of chosen item
      R[i] <- A[chosen_idx]

      # Remove chosen
      A <- A[-chosen_idx]
      Gamma <- Gamma[-chosen_idx]

      # Renormalize
      s <- sum(Gamma)
      if (s > 0) {
        Gamma <- Gamma / s
      } else {
        # If all remaining mass is 0, pick deterministically in order
        Gamma <- c(1, rep(0, length(Gamma) - 1L))
      }
    }

    # Last remaining item
    R[t] <- A
    rank_matrix[j, ] <- R
  }

  out <- as.data.frame(rank_matrix, stringsAsFactors = FALSE)
  colnames(out) <- ordinal_seq(t)
  out
}
