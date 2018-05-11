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
    group_by(term) %>%
    mutate(odds_ratio = exp(estimate),
           ciub = exp(estimate + 2 * std.error.naive),
           cilb = exp(estimate - 2 * std.error.naive)) %>%
    summarize(ub = max(odds_ratio),
              lb = min(odds_ratio),
              ciub = max(ciub),
              cilb = min(cilb))

  base_pd <- x$base_case %>%
    filter(controls %in% rad_ctls) %>%
    mutate(odds_ratio = exp(estimate),
           base_ciub = exp(estimate + 2 * std.error.naive),
           base_cilb = exp(estimate - 2 * std.error.naive))

  pd <- merge(base_pd, sens_pd, by = "term")

  if (include_benchmark) {
    if ("bm" %in% x$base_case$method) {
      bm_pd <- x$base_case[!(x$base_case$controls %in% rad_ctls), ]
      bm_pd <- bm_pd %>%
        mutate(odds_ratio = exp(estimate),
               base_ciub = exp(estimate + 2 * std.error),
               base_cilb = exp(estimate - 2 * std.error))

      pd <- bind_rows(bm_pd, pd)

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
                alpha = .4, na.rm = include_benchmark) +
    geom_linerange(aes(ymin = lb, ymax = ub), na.rm = include_benchmark,
                   size = 4) +
    geom_errorbar(aes(ymin = base_cilb, ymax = base_ciub),
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


#' Plot risk-treatment and calibration plots for a \code{policy}
#'
#' @param x object of class \code{policy}
#' @param nbins number of bins to use (e.g., for calibration plots)
#' @param method passed to \code{\link[ggplot2]{geom_smooth}}
#' @param method.args passed to \code{\link[ggplot2]{geom_smooth}}
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return a named list of ggplot object
#' @export
plot.policy <- function(x,
                        nbins = 10,
                        method = "glm",
                        method.args = list(family = "quasibinomial"),
                        ...) {
  s_treatment <- x$treatment
  v_outcome <- rlang::sym(x$outcome)
  v_group <- rlang::sym(x$grouping)
  v_treatment <- rlang::sym(s_treatment)
  v_riskcol <- rlang::sym(x$risk_col)

  # Risk vs. treatment
  p_risk_v_trt <-
    ggplot(data = x$data, aes_string(
      x = paste0(x$risk_col, "__"),
      y = paste0("as.numeric(", x$treatment, ")")
    )) +
    geom_smooth(aes_string(color = x$grouping),
                method = method, method.args = method.args) +
    scale_x_continuous(paste0("\nEstimated risk (", x$risk_col, ")"),
                       labels = scales::percent) +
    scale_y_continuous(paste0("Outcome (", x$outcome, ")"), labels = scales::percent)

  # Risk calibration plots, conditioned on control/treatment
  calib_ctl_pd <- x$data %>%
    filter(as.numeric(!!v_treatment) == 0, fold__ == "test") %>%
    mutate(risk_bin = ntile(resp_ctl__, 10)) %>%
    group_by(risk_bin, !!v_group) %>%
    summarize(N = n(),
              pout = mean(!!v_outcome),
              mresp = mean(resp_ctl__)) %>%
    mutate(type = sprintf("Outcome given %s = 0", s_treatment))

  calib_trt_pd <- x$data %>%
    filter(as.numeric(!!v_treatment) == 1, fold__ == "test") %>%
    mutate(risk_bin = ntile(resp_trt__, 10)) %>%
    group_by(risk_bin, !!v_group) %>%
    summarize(N = n(),
              pout = mean(!!v_outcome),
              mresp = mean(resp_trt__)) %>%
    mutate(type = sprintf("Outcome given %s = 1", s_treatment))

  calib_pd <- bind_rows(calib_ctl_pd, calib_trt_pd)

  p_calib <- ggplot(calib_pd, aes(x = mresp, y = pout)) +
    geom_abline(intercept = 0, slope = 1, size = 2,
                color = theme_get()$panel.grid.minor$colour) +
    geom_point(aes_string(size = "N", color = x$grouping)) +
    scale_size_area() +
    scale_x_continuous(sprintf("\nEstimated risk (%s)", x$risk_col),
                       labels = scales::percent) +
    scale_y_continuous(sprintf("Proportion of %s = 1\n", x$outcome),
                       labels = scales::percent) +
    facet_grid(. ~ type) +
    theme(legend.position = "top")

  ret <- list(
    risk_v_trt = p_risk_v_trt,
    calibration = p_calib
  )

  ret
}


#' Plot changes in risk-treatment for a \code{sensitive_policy}
#'
#' @param x object of class \code{sensitive_policy}
#' @param down_sample (Optional) proportion (between 0 and 1) or number (greater
#'   than 1) of rows to sample from each group, if down sampling the data;
#'   default is 30
#' @param method passed to \code{\link[ggplot2]{geom_smooth}}
#' @param method.args passed to \code{\link[ggplot2]{geom_smooth}}
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return a ggplot object
#' @export
plot.sensitive_policy <- function(x,
                                  down_sample = 30,
                                  method = "glm",
                                  method.args = list(family = "quasibinomial"),
                                  ...) {
  s_treatment <- x$treatment
  v_outcome <- rlang::sym(x$outcome)
  v_group <- rlang::sym(x$grouping)
  v_treatment <- rlang::sym(s_treatment)
  v_riskcol <- rlang::sym(x$risk_col)

  vanilla_df <- .down_sample(x$data %>% group_by(!!v_group),
                             down_sample,
                             verbose = FALSE) %>%
    mutate(weights__ = 1, type = "Original", !!v_treatment := ptrt__)

  caption <- paste0("Points are down-sampled to ", nrow(vanilla_df), "/",
                    nrow(x$data), " rows (",
                    format(nrow(vanilla_df)/nrow(x$data) * 100), "%)")

  sampled_ids <- vanilla_df$id_sens__

  sens_df <- x$sens_data %>%
    filter(id_sens__ %in% sampled_ids) %>%
    mutate(type = "Sensitized")

  pd <- rbind(x$data %>%
                select(!!x$treatment, !!paste0(x$risk_col, "__"),
                       !!x$grouping, risk__, ptrt__) %>%
                mutate(weights__ = 1, type = "Original"),
              x$sens_data %>%
                select(!!x$treatment, !!paste0(x$risk_col, "__"),
                              !!x$grouping, risk__, ptrt__, weights__) %>%
                mutate(type = "Sensitized"))

  # Risk vs. treatment
  ret <-
    ggplot(data = vanilla_df, aes_string(
      x = "risk__",
      # x = paste0(x$risk_col, "__"),
      y = x$treatment,
      alpha = "type",
      color = x$grouping
    )) +
    geom_smooth(formula = y ~ x, data = pd, se = FALSE, method = method,
                method.args = method.args,
                aes(weight = weights__, linetype = type)) +
    geom_point(aes(shape = type, size = weights__)) +
    geom_point(data = sens_df, aes(shape = type, size = weights__)) +
    geom_line(data = sens_df, aes(group = id_sens__),
              size = .2, linetype = "dashed", alpha = .3) +
    scale_size_area(max_size = 3) +
    scale_linetype_manual("Type",
                          limits = c("Original", "Sensitized"),
                          values = c("dotted", "solid")) +
    scale_shape_manual("Type",
                       limits = c("Original", "Sensitized"),
                       values = c(1, 19)) +
    scale_alpha_manual("Type",
                       limits = c("Original", "Sensitized"),
                       values = c(.4, .8)) +
    scale_x_continuous(paste0("\nEstimated risk (", x$risk_col, ")")) +
    scale_y_continuous("Estimated probability of treatment\n",
                       labels = scales::percent) +
    labs(caption = caption)

  ret
}
