plot.optimsens <- function(os, controls = NULL) {
  base_pd <- os$base_case
  base_pd$odds_ratio <- exp(base_pd$estimate)
  sens_pd <- os$results %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(ciub = estimate + 2 * std.error.naive,
                  cilb = estimate - 2 * std.error.naive) %>%
    dplyr::summarize(
      ub = exp(max(estimate)),
      lb = exp(min(estimate)),
      ciub = exp(max(ciub)),
      cilb = exp(min(cilb))
    )

  pd <- merge(base_pd, sens_pd, by = "term")

  ggplot(pd, aes(x = term, y = odds_ratio)) +
    geom_errorbar(aes(ymin = cilb, ymax = ciub), size = .5, width = .2,
                alpha = .4) +
    geom_linerange(aes(ymin = lb, ymax = ub), size = 4) +
    geom_errorbar(aes(ymin = cilb, ymax = ciub),
                  size = .8, width = 0, alpha = .6) +
    geom_point(color = "black", fill = "white", size = 3, shape = 21) +
    scale_y_continuous(sprintf("Odds of treatment (v. %s)\n", os$base_group)) +
    scale_x_discrete(element_blank())
}
