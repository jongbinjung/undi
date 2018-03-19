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
#' @param down_sample (Optional) proportion (between 0 and 1) or number (greater
#'   than 1) of rows to sample, if down sampling the (test) data; default is 1
#'   (i.e., use all data)
#' @param seed random seed to set
#'
#' @return tidy data frame of rad coefficients
#'
#' @export
compute_rad <-
  function(pol,
           controls = NULL,
           base_group = NULL,
           minority_groups = NULL,
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

    if (is.null(base_group)) {
      if (is.factor(group_col)) {
        base_group <- levels(group_col)[1]
      } else {
        base_group <- unique(group_col)[1]
      }
    }

    if (is.null(minority_groups)) {
      if (is.factor(group_col)) {
        minority_groups <- levels(group_col)[-1]
      } else {
        minority_groups <- unique(group_col)[-1]
      }
    }

    test_df <- .down_sample(d[d$fold__ == "test", ], down_sample)

    ret <- purrr::map_dfr(minority_groups, function(comp) {
      target_group_ind <- test_df[[pol$grouping]] %in% c(base_group, comp)

      tmp_df <- test_df[target_group_ind, ]
      tmp_df[[pol$grouping]] <- forcats::fct_drop(tmp_df[[pol$grouping]])
      tmp_df[[pol$grouping]] <- forcats::fct_relevel(tmp_df[[pol$grouping]],
                                                   base_group)

      coefs <- .pull_coefs(tmp_df, pol$treatment, pol$grouping,
                           c("risk__", controls),
                           fun = pol$fit2)

      coefs[grepl(pol$grouping, coefs$term), ]
    })

    return(ret)
  }

#' Compute nonparametric risk-adjusted disparate impact estimate for a policy
#'
#' @param pol object of class policy
#' @param controls character vector of additional controls to consider in the
#'   second-stage model
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#'
#' @details Current implementation uses \code{loess} as a "non-parametric"
#'   estimate
#'
#' @return data frame of non-parametric estimates of average treatment for each
#'   group, assuming the entire population is one of each group
#'
#' @export
compute_nprad <-
  function(pol,
           controls = NULL,
           base_group = NULL,
           minority_groups = NULL) {
    # Input validation
    if (!("policy" %in% class(pol))) {
      stop("Expected object of class policy")
    }

    if (length(base_group) > 1) {
      stop("Specify a single base group.\n\tGot: ", base_group)
    }

    if (!is.null(controls)) {
      stop("controls are not yet implemented for non-parametric rad")
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

    if (is.null(base_group)) {
      if (is.factor(group_col)) {
        base_group <- levels(group_col)[1]
      } else {
        base_group <- unique(group_col)[1]
      }
    }

    if (is.null(minority_groups)) {
      if (is.factor(group_col)) {
        minority_groups <- levels(group_col)[-1]
      } else {
        minority_groups <- unique(group_col)[-1]
      }
    }

    test_df <- d[d$fold__ == "test", ]

    # Fit loess models for each group (in parallel)
    loess_ms <-
      foreach::foreach(group = c(base_group, minority_groups)) %dopar% {
        group_ind <- test_df[[pol$grouping]] == group

        train_df <- test_df[group_ind, ]
        # Fit on logit scale to avoid values out-of-bound
        list(model = stats::loess(logit(ptrt__) ~ risk__, data = train_df),
             group = group)
    }

    # Use loess model to predict average treatment, assuming the entire
    # population is from a single group
    ret <- purrr::map_dfr(loess_ms, function(loess_m) {
      model <- loess_m[["model"]]
      group <- loess_m[["group"]]

      tmp_df <- test_df
      tmp_df[[pol$grouping]] <- group

      dplyr::tibble(!!pol$grouping := group,
                    ptrt =  mean(inv_logit(stats::predict(model, data = tmp_df))))
    })

    # Compute odds and odds ratio v. base_group
    ret$odds <- ret$ptrt / (1 - ret$ptrt)
    base_odds <- ret[ret[[pol$grouping]] == base_group, ]$odds

    ret$or <- ret$odds / base_odds

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
    members <- unique(group_col)

    check_groups <- sapply(c(base_group, minority_groups),
                           function(x) x %in% members)
    if (!all(check_groups)) {
      stop(sprintf("%s - not members of %s",
                   paste0(c(base_group, minority_groups)[!check_groups],
                          collapse = ","),
                   pol$grouping))
    }

    if (is.null(base_group)) {
      if (is.factor(group_col)) {
        base_group <- levels(group_col)[1]
      } else {
        base_group <- unique(group_col)[1]
      }
    }

    if (is.null(minority_groups)) {
      if (is.factor(group_col)) {
        minority_groups <- levels(group_col)[-1]
      } else {
        minority_groups <- unique(group_col)[-1]
      }
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

      coefs <- .pull_coefs(tmp_df, pol$treatment, pol$grouping,
                           c(controls),
                           fun = function(f, d, w)
                             stats::glm(f, d, family = stats::binomial))

      coefs[grepl(pol$grouping, coefs$term), ]
    })

    if (kitchen_sink) ret$controls <- "kitchen sink"

    return(ret)
  }
