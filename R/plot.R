#' Plot min/max sensitivity results for a policy from a \code{sens} object
#'
#' @param x object of class \code{sens}
#' @param include_benchmark logical; whether to include the two extreme
#'   benchmark test results; if TRUE, the returned ggplot object is facetted by
#'   test type (default = TRUE)
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return a ggplot object
#' @export
plot.sens <- function(x, include_benchmark = TRUE, ...) {
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
    if ("bm" %in% x$base_case$method) {
      bm_pd <- x$base_case[!(x$base_case$controls %in% rad_ctls), ]
      bm_pd <- bm_pd %>%
        dplyr::mutate(ciub = exp(estimate + 2 * std.error),
                      cilb = exp(estimate - 2 * std.error))
      bm_pd$odds_ratio <- exp(bm_pd$estimate)

      pd <- dplyr::bind_rows(bm_pd, pd)

      pd$controls <- forcats::fct_inorder(pd$controls)
    } else {
      warning(paste("include_benchmark = TRUE;",
                    "but benchmark results not computed in", class(x)[1],
                    "object\n",
                    "rerun", class(x)[1], "with include_benchmark = TRUE if",
                    "you want to plot benchmark results"))
      include_benchmark <- FALSE
    }
  }

  p <- ggplot(pd, aes(x = term, y = odds_ratio)) +
    geom_hline(yintercept = 1, size = 2,
               color = ggplot2::theme_get()$panel.grid.minor$colour) +
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


#' Plot min/max sensitivity results for a policy from a \code{sens} object
#'
#' @param x object of class \code{policy}
#' @param nbins number of bins to use (e.g., for calibration plots)
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return a named list of ggplot object
#' @export
plot.policy <- function(x, nbins = 10, ...) {
  v_outcome <- rlang::sym(x$outcome)
  v_group <- rlang::sym(x$grouping)
  v_treatment <- rlang::sym(x$treatment)
  v_riskcol <- rlang::sym(x$risk_col)

  # Risk vs. treatment
  p_risk_v_trt <-
    ggplot(data = x$data, aes_string(
      x = paste0(x$risk_col, "__"),
      y = paste0("as.numeric(", x$treatment, ")")
    )) +
    geom_smooth(aes_string(color = x$grouping),
                method = "glm", method.args = list(family = "binomial")) +
    scale_x_continuous(paste0("\nEstimated risk (", x$risk_col, ")"),
                       labels = scales::percent) +
    scale_y_continuous(paste0("Outcome (", x$outcome, ")"), labels = scales::percent)

  # Risk calibration plots, conditioned on control/treatment
  calib_ctl_pd <- x$data %>%
    dplyr::filter(!!v_treatment == 0) %>%
    dplyr::mutate(risk_bin = dplyr::ntile(resp_ctl__, 10)) %>%
    dplyr::group_by(risk_bin, !!v_group) %>%
    dplyr::summarize(N = n(),
                     pout = mean(!!v_outcome),
                     mresp = mean(resp_ctl__)) %>%
    dplyr::mutate(type = "Outcome given control")

  calib_trt_pd <- x$data %>%
    dplyr::filter(!!v_treatment == 1) %>%
    dplyr::mutate(risk_bin = dplyr::ntile(resp_trt__, 10)) %>%
    dplyr::group_by(risk_bin, !!v_group) %>%
    dplyr::summarize(N = n(),
                     pout = mean(!!v_outcome),
                     mresp = mean(resp_trt__)) %>%
    dplyr::mutate(type = "Outcome given treatment")

  calib_pd <- dplyr::bind_rows(calib_ctl_pd, calib_trt_pd)

  p_calib <- ggplot(calib_pd, aes(x = mresp, y = pout)) +
    geom_abline(intercept = 0, slope = 1, size = 2,
                color = theme_get()$panel.grid.minor$colour) +
    geom_point(aes(size = N, color = x)) +
    scale_size_area() +
    scale_x_continuous("\nEstimated risk", labels = scales::percent) +
    scale_y_continuous("Proportion of outcome = 1\n", labels = scales::percent) +
    facet_grid(. ~ type) +
    theme(legend.position = "top")

  ret <- list(
    risk_v_trt = p_risk_v_trt,
    calibration = p_calib
  )
}
