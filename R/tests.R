#' Compute risk-adjusted disparate impact estimate for a policy
#'
#' @param controls character vector of additional controls to consider in the
#'   second-stage model
#' @param down_sample (Optional) proportion (between 0 and 1) or number (greater
#'   than 1) of rows to sample, if down sampling the (test) data; default is 1
#'   (i.e., use all data)
#' @param seed random seed to set
#' @param ... additional arguments to pass to \code{fit} function from
#'   \code{\link{rad_control}} for fine-tuning
#'
#' @return tidy data frame with columns \item{term}{the group members considered
#'   minority} \item{estimate}{log-odds of treatment, relative to base_group
#'   (equivalent to logistic regression coefficient)}
#'   \item{std.error}{coefficient standard errors for \code{*_coef} methods; TO
#'   BE IMPLEMENTED for \code{*_avg} methods (FOR NOW, ALL ZERO)!}
#'   \item{statistic/p.value}{(for \code{*_coef} methods) corresponding values
#'   from model fit} \item{ptrt_base/minor}{(for \code{*_avg} methods) estimated
#'   average treatment probability for base/minority groups} \item{method}{the
#'   method used} \item{controls}{features controlled for}
#'
#' @inheritParams .validate_input
#' @inheritParams rad_control
#' @export
compute_rad <-
  function(pol,
           controls = NULL,
           base_group = NULL,
           minority_groups = NULL,
           fit_fn = "logit_coef",
           down_sample = 1,
           use_speedglm = TRUE,
           seed = round(stats::runif(1)*1e4),
           ...) {
    set.seed(seed)
    # Input validation
    groups <- .validate_input(pol, base_group, minority_groups)
    base_group <- groups$base
    minority_groups <- groups$minority

    rc <- rad_control(pol, controls, fit_fn = fit_fn,
                      use_speedglm = use_speedglm)

    d <- pol$data

    # Restrict data to groups of interest
    target_group_ind <- d[[pol$grouping]] %in% c(base_group, minority_groups)
    d <- d[target_group_ind, ]
    d[[pol$grouping]] <- forcats::fct_drop(d[[pol$grouping]])

    # Down-sample and filter to test fold
    test_df <- .down_sample(d[d$fold__ == "test", ], down_sample)

    # Make sure that the base_group is first level
    # (see @details of .compute_estimate)
    test_df[[pol$grouping]] <- forcats::fct_relevel(test_df[[pol$grouping]],
                                                    base_group)

    ret <- .compute_estimate(test_df, rc, ...)
    return(ret)
  }


#' Compute benchmark test for disparate impact of a policy
#'
#' @param controls character vector of additional controls to consider (i.e.,
#'   valid benchmarks)
#' @param kitchen_sink logical; if TRUE, ignore \code{controls} argument, and
#'   include all variables given with the policy, i.e., policy$features
#'   (default: FALSE)
#'
#' @return tidy data frame of benchmark coefficients
#'
#' @inheritParams .validate_input
#' @export
compute_bm <-
  function(pol,
           controls = NULL,
           kitchen_sink = FALSE,
           base_group = NULL,
           minority_groups = NULL) {
    # Input validation
    groups <- .validate_input(pol, base_group, minority_groups)
    base_group <- groups$base
    minority_groups <- groups$minority

    if (kitchen_sink) {
      controls <- pol$features
    }

    d <- pol$data
    test_df <- d[d$fold__ == "test", ]

    ret <- purrr::map_dfr(minority_groups, function(comp) {
      target_group_ind <- test_df[[pol$grouping]] %in% c(base_group, comp)

      tmp_df <- test_df[target_group_ind, ]
      tmp_df[[pol$grouping]] <- forcats::fct_drop(tmp_df[[pol$grouping]])
      tmp_df[[pol$grouping]] <- forcats::fct_relevel(tmp_df[[pol$grouping]],
                                                   base_group)

      coefs <- .get_estimate(tmp_df, pol$treatment, pol$grouping,
                             c(controls),
                             fun = function(f, d, w)
                               stats::glm(f, d, family = stats::binomial))

      coefs[grepl(pol$grouping, coefs$term), ]
    })

    if (kitchen_sink) ret$controls <- "kitchen sink"

    return(ret)
  }



#' Compute outcome tests for disparate impact of a policy
#'
#' @param controls character vector of additional controls to consider (i.e.,
#'   conditional groupings)
#'
#' @return data frame, grouped by the group and additional columns specified in
#'   the \code{controls} arguments and corresponding \code{hitrate}, i.e.,
#'   \code{P(outcome = 1 | treatment = risk_treatment)}, where
#'   \code{risk_treatment = ifelse(pol$risk_col == "resp_trt", 1, 0)}.
#'
#' @inheritParams .validate_input
#' @export
compute_ot <- function(pol,
                       controls = NULL) {
  d <- pol$data
  test_df <- d[d$fold__ == "test", ]

  v_treatment <- rlang::sym(pol$treatment)
  v_outcome <- rlang::sym(pol$outcome)
  risk_treatment <- ifelse(pol$risk_col == "resp_trt", 1, 0)

  test_df %>%
    filter(!!v_treatment == risk_treatment) %>%
    group_by_(.dots = c(pol$grouping, controls)) %>%
    summarize(hitrate = mean(!!v_outcome), count = n())
}

#' Given a data frame and \code{\link{rad_control} object}, compute RAD estimate
#'
#' @param d data frame that has all necessary columns
#' @param rc a \code{\link{rad_control}} object
#' @param ... additional arguments passed to \code{rc$fit} function
#'
#' @details This helper method relies on the proper factoring of data \code{d}.
#'   The grouping variable is extracted from \code{rc$grouping}, while the
#'   \code{base_group} and \code{minority_groups} are determined by the levels
#'   ordering in \code{d[[rc$grouping]]}, i.e., make sure groups are properly
#'   ordered before calling!
#'
#' @return tidy dataframe of estimated rad results
.compute_estimate <- function(d,  rc, weighted = FALSE, ...) {
  # Fit model
  groups <- .get_groups(d[[rc$grouping]])
  base_group <- groups[1]
  minority_groups <- groups[-1]

  ret <- purrr::map_dfr(minority_groups, function(group) {
    target_group_ind <- d[[rc$grouping]] %in% c(base_group, group)

    tmp_df <- d[target_group_ind, ]

    if (weighted) {
      m <- rc$fit(tmp_df, w = tmp_df$weights__, ...)
    } else {
      m <- rc$fit(tmp_df, ...)
    }

    if (rc$method == "coef") {
      coefs <- m %>%
        broom::tidy() %>%
        mutate(method = rc$method, controls = rc$label)

      coefs[grepl(rc$grouping, coefs$term), ]
    } else if (rc$method == "avg") {
      ptrt <- purrr::map_dbl(c(base_group, group), function(x) {
        counter_df <- tmp_df
        counter_df[[rc$grouping]] <- x
        if (weighted) {
          stats::weighted.mean(rc$pred(m, counter_df), w = counter_df$weights__)
        } else {
          mean(rc$pred(m, counter_df))
        }
      })

      odds <- ptrt / (1 - ptrt)
      or <- odds[2]/odds[1]

      # TODO: estimate standard errors? (just create column of 0 for now)
      tibble(term = paste0(rc$grouping, group),
                    estimate = log(or),
                    std.error = 0,
                    ptrt_base = ptrt[1],
                    ptrt_minor = ptrt[2],
                    method = rc$method,
                    controls = rc$label)
    } else {
      stop("Unknown method specification from rad_control:", rc$method)
    }
  })
  }


#' Validate input for
#'
#' @param pol object of class policy
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#'
#' @return list of validated group members \item{base}{base group
#'   members}\item{minority}{minority group members}
.validate_input <- function(pol, base_group, minority_groups) {
  if (!("policy" %in% class(pol))) {
    stop("Expected object of class policy")
  }

  if (length(base_group) > 1) {
    stop("Specify a single base group.\n\tGot: ", base_group)
  }

  d <- pol$data

  group_col <- d[[pol$grouping]]
  groups <- .get_groups(group_col)

  check_groups <- sapply(c(base_group, minority_groups),
                         function(x) x %in% groups)
  if (!all(check_groups)) {
    stop(sprintf("%s - not member of %s",
                 paste0(c(base_group, minority_groups)[!check_groups],
                        collapse = ","),
                 pol$grouping))
  }

  if (is.null(base_group)) {
    base_group <- groups[1]
  }

  if (is.null(minority_groups)) {
    minority_groups <- groups[!(groups == base_group)]
  }

  list(base = base_group, minority = minority_groups)
}
