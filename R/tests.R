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
#'
#' @return tidy data frame of rad coefficients
#'
#' @export
compute_rad <-
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
                             stats::glm(f, d, family = binomial))

      coefs[grepl(pol$grouping, coefs$term), ]
    })

    if (kitchen_sink) ret$controls <- "kitchen sink"

    return(ret)
  }
