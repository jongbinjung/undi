#' Plot optimized sensitivity results for a policy
#'
#' @param x object of class \code{optimsens}
#' @param include_benchmark logical; whether to include the two extreme
#'   benchmark test results; if TRUE, the returned ggplot object is facetted by
#'   test type (default = TRUE)
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return a ggplot object
#' @export
plot.optimsens <- function(x, include_benchmark = TRUE, ...) {
  rad_ctls <- unique(x$results$controls)
  sens_pd <- x$results %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(ciub = estimate + 2 * std.error.naive,
                  cilb = estimate - 2 * std.error.naive) %>%
    dplyr::summarize(
      ub = exp(max(estimate)),
      lb = exp(min(estimate)),
      ciub = exp(max(ciub)),
      cilb = exp(min(cilb))
    )

  base_pd <- x$base_case[x$base_case$controls %in% rad_ctls, ]
  base_pd$odds_ratio <- exp(base_pd$estimate)

  pd <- merge(base_pd, sens_pd, by = "term")

  if (include_benchmark) {
    bm_pd <- x$base_case[!(x$base_case$controls %in% rad_ctls), ]
    bm_pd <- bm_pd %>%
      dplyr::mutate(ciub = exp(estimate + 2 * std.error),
                    cilb = exp(estimate - 2 * std.error))
    bm_pd$odds_ratio <- exp(bm_pd$estimate)

    pd <- dplyr::bind_rows(bm_pd, pd)

    pd$controls <- forcats::fct_inorder(pd$controls)
  }

  p <- ggplot(pd, aes(x = term, y = odds_ratio)) +
    geom_errorbar(aes(ymin = cilb, ymax = ciub), size = .5, width = .2,
                alpha = .4) +
    geom_linerange(aes(ymin = lb, ymax = ub), size = 4) +
    geom_errorbar(aes(ymin = cilb, ymax = ciub),
                  size = .8, width = 0, alpha = .6) +
    geom_point(color = "black", fill = "white", size = 3, shape = 21) +
    scale_y_continuous(sprintf("Odds of treatment (v. %s)\n", x$base_group)) +
    scale_x_discrete(element_blank())

  if (include_benchmark) {
    p <- p +
      facet_grid(~ method + controls)
  }

  return(p)
}
