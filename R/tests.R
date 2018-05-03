#' Compute risk-adjusted disparate impact estimate for a policy
#'
#' @param pol object of class policy
#' @param controls character vector of additional controls to consider in the
#'   second-stage model
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#' @param fit_fn string indicating the fitting proceedure used.
#' @param down_sample (Optional) proportion (between 0 and 1) or number (greater
#'   than 1) of rows to sample, if down sampling the (test) data; default is 1
#'   (i.e., use all data)
#' @param seed random seed to set
#' @param ... additional arguments to pass to \code{fit} function from
#'   \code{\link{di_model}} for fine-tuning
#'
#' @return tidy data frame of rad coefficients
#'
#' @export
compute_rad <-
  function(pol,
           controls = NULL,
           base_group = NULL,
           minority_groups = NULL,
           fit_fn = c("logit", "gam"),
           down_sample = 1,
           seed = round(stats::runif(1)*1e4),
           ...) {
    set.seed(seed)
    # Input validation
    if (!("policy" %in% class(pol))) {
      stop("Expected object of class policy")
    }

    if (length(base_group) > 1) {
      stop("Specify a single base group.\n\tGot: ", base_group)
    }

    fit_fn <- match.arg(fit_fn)

    dm <- di_model(pol, controls, fit_fn = fit_fn)

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
      minority_groups <- groups[-1]
    }

    # Restrict data to groups of interest
    target_group_ind <- d[[pol$grouping]] %in% c(base_group, minority_groups)
    d <- d[target_group_ind, ]
    d[[pol$grouping]] <- forcats::fct_drop(d[[pol$grouping]])

    # Down-sample and filter to test fold
    test_df <- .down_sample(d[d$fold__ == "test", ], down_sample)

    # Make sure that the base_group is first level
    test_df[[pol$grouping]] <- forcats::fct_relevel(test_df[[pol$grouping]],
                                                    base_group)

    ret <- .compute_estimate(test_df, pol$grouping, dm, ...)

    return(ret)
  }


#' Compute benchmark test for disparate impact of a policy
#'
#' @param pol object of class policy
#' @param controls character vector of additional controls to consider (i.e.,
#'   valid benchmarks)
#' @param kitchen_sink logical; if TRUE, ignore \code{controls} argument, and
#'   include all variables given with the policy, i.e., policy$features
#'   (default: FALSE)
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#'
#' @return tidy data frame of benchmark coefficients
#'
#' @export
compute_bm <-
  function(pol,
           controls = NULL,
           kitchen_sink = FALSE,
           base_group = NULL,
           minority_groups = NULL) {
    # Input validation
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
      minority_groups <- groups[-1]
    }

    if (kitchen_sink) {
      controls <- pol$features
    }

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


#' Compute risk-adjusted disparate impact estimate for a policy
#'
#' @param pol object of class policy
#' @param controls character vector of additional controls to consider in the
#'   second-stage model
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#' @param difun a "disparate impact function" of the form \code{f(data, weights
#'   = NULL)} used for predicting treatment conditional on group membership;
#'   using \code{glm} with \code{family = quasibinomial} by default; the
#'   \code{weights} argument is only used for sensitivity and *must* be
#'   initialized to \code{NULL} (or the equivalent of non-weighted fitting)
#' @param down_sample (Optional) proportion (between 0 and 1) or number (greater
#'   than 1) of rows to sample, if down sampling the (test) data; default is 1
#'   (i.e., use all data)
#' @param seed random seed to set
#'
#' @return tidy data frame of rad coefficients
#'
#' @export
compute_rad_old <-
  function(pol,
           controls = NULL,
           base_group = NULL,
           minority_groups = NULL,
           difun = NULL,
           down_sample = 1,
           seed = round(stats::runif(1)*1e4)) {
    set.seed(seed)
    # Input validation
    if (!("policy" %in% class(pol))) {
      stop("Expected object of class policy")
    }

    if (length(base_group) > 1) {
      stop("Specify a single base group.\n\tGot: ", base_group)
    }

    if (is.null(difun)) {
      difun <- models()$glm$fit
    }

    d <- pol$data

    group_col <- d[[pol$grouping]]
    members <- unique(group_col)

    check_groups <- sapply(c(base_group, minority_groups),
                           function(x) x %in% members)
    if (!all(check_groups)) {
      stop(sprintf("%s - not members of %s",
                   paste0(c(base_group, minority_groups)[!check_groups],
                          collapse = ","),
                   pol$grouping))
    }

    if (is.null(base_group) | is.null(minority_groups)) {
      groups <- .get_groups(group_col)
    }

    if (is.null(base_group)) {
      base_group <- groups[1]
    }

    if (is.null(minority_groups)) {
      minority_groups <- groups[-1]
    }

    test_df <- .down_sample(d[d$fold__ == "test", ], down_sample)

    ret <- purrr::map_dfr(minority_groups, function(comp) {
      target_group_ind <- test_df[[pol$grouping]] %in% c(base_group, comp)

      tmp_df <- test_df[target_group_ind, ]
      tmp_df[[pol$grouping]] <- forcats::fct_drop(tmp_df[[pol$grouping]])
      tmp_df[[pol$grouping]] <- forcats::fct_relevel(tmp_df[[pol$grouping]],
                                                   base_group)

      coefs <- .get_estimate(tmp_df, pol$treatment, pol$grouping,
                             c("risk__", controls),
                             fun = difun)

      coefs[grepl(pol$grouping, coefs$term), ]
    })

    return(ret)
  }


#' Fit second stage model and compute disparate impact estimate
#'
#' Given a data frame, group column name, and functions for fitting and
#' predicting treatment, report the average odds ratio of treatment between
#' groups.
#'
#' @param d data frame that has all necessary columns
#' @param cn_group character string of grouping column name
#' @param dm a \code{\link{di_model}} object
#' @param alt_fit alternative \code{fit} function to use; \code{dm$fit} is used
#'   if \code{NULL}; general purpose is to allow for weights
#' @param ... additional arguments passed to \code{fit} function
#'
#' @return tidy dataframe of the glm model
.compute_estimate <-
  function(d,
           cn_group,
           dm,
           weighted = FALSE,
           ...) {
  # Fit model
  groups <- .get_groups(d[[cn_group]])
  base_group <- groups[1]
  minority_groups <- groups[-1]

  ret <- purrr::map_dfr(minority_groups, function(group) {
    target_group_ind <- d[[cn_group]] %in% c(base_group, group)

    tmp_df <- d[target_group_ind, ]

    if (weighted) {
      m <- dm$fit(tmp_df, w = tmp_df$weights__, ...)
    } else {
      m <- dm$fit(tmp_df, ...)
    }

    ptrt <- purrr::map_dbl(c(base_group, group), function(x) {
      counter_df <- tmp_df
      counter_df[[cn_group]] <- x
      if (weighted) {
        stats::weighted.mean(dm$pred(m, counter_df), w = counter_df$weights__)
      } else {
        mean(dm$pred(m, counter_df))
      }
    })

    odds <- ptrt / (1 - ptrt)
    or <- odds[2]/odds[1]

    # TODO: estimate standard errors? (just create column of 0 for now)
    dplyr::tibble(term = paste0(cn_group, group),
                  estimate = or,
                  std.error = 0,
                  ptrt_base = ptrt[1],
                  ptrt_minor = ptrt[2],
                  controls = dm$label)
  })
}
