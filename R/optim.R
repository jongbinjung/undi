#' Find parameter values that min/maximize undi results
#'
#' Within specified range of sensitivity parameters, find the ones that achieve
#' minimum/maximum undi results
#'
#' @param u object of class \code{undi}
#' @param range_q 2D vector specifying min/max value of p(u = 1 | x)
#' @param range_dp 2D vector specifying min/max value of change in log-odds of
#'   treat = 1 if u = 1
#' @param range_d0 2D vector specifying min/max value of change in log-odds of
#'   response = 1 if treat = 0 and u = 1
#' @param range_d1 2D vector specifying min/max value of change in log-odds of
#'   response = 1 of treat = 1 and u = 1
#' @param (Optional) base_group single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param (Optional) minority_groups groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#' @param debug logical flag, if TRUE, returns a list of results and the
#'   expanded data frame used to fit model
#'
#' @details All sensitivity parameters (\code{q, dp, d0, d1}) can be provided in
#'   one of three formats, determined by the \code{length} of each argument:
#'   \describe{ \item{if \code{length(arg) = 1}}{single value applied to all
#'   observations (rows)} \item{if \code{length(arg) = }number of levels in
#'   grouping variable}{each parameter setting applied to corresponding level in
#'   group} \item{if \code{length(arg) = nrow(u$data)}}{each parameter applied
#'   to corresponding rows}} Note that if \code{compare} is specified, the
#'   number of grouping levels is effectively the length of \code{compare}
#'
#' @return \code{tidy} dataframe of second-stage model coefficients after
#'   applying sensitivity parameters, with a nested column of sensitivity params
#'
#' @export
optimsens <-
  function(u,
           range_q = c(0, 1),
           range_dp = c(0, log(2)),
           range_d0 = c(0, log(2)),
           range_d1 = c(0, log(2)),
           optim_params = NULL,
           base_group = NULL,
           minority_groups = NULL,
           debug = FALSE) {
  # Input validation
  if (!("undi" %in% class(u))) {
    stop("Expected object of class undi")
  }

  if (length(base_group) > 1) {
    stop("Specify a single base group.\n\tGot: ", base_group)
  }

  # range_q = c(0, 1)
  # range_dp = c(0, log(2))
  # range_d0 = c(0, log(2))
  # range_d1 = c(0, log(2))

  group_col <- u$data[[u$grouping]]
  members <- unique(group_col)

  check_groups <- sapply(c(base_group, minority_groups),
                         function(x) x %in% members)
  if (!all(check_groups)) {
    stop(sprintf("%s - not members of %s",
                 paste0(c(base_group, minority_groups)[!check_groups],
                        collapse = ","),
                 u$grouping))
  }

  # Optimization hyper parameters
  params_lower <- c(
    range_q[1],   # Min P(u = 1 | base )
    range_q[1],   # Min P(u = 1 | minority )
    range_dp[1],  # Min alpha_base
    range_dp[1],  # Min alpha_minority
    range_d0[1],  # Min delta0_base
    range_d0[1],  # Min delta0_minority
    range_d1[1],  # Min delta1_base
    range_d1[1]   # Min delta1_minority
  )
  params_upper <- c(
    range_q[2],   # Max P(u = 1 | base )
    range_q[2],   # Max P(u = 1 | minority )
    range_dp[2],  # Max alpha_base
    range_dp[2],  # Max alpha_minority
    range_d0[2],  # Max delta0_base
    range_d0[2],  # Max delta0_minority
    range_d1[2],  # Max delta1_base
    range_d1[2]   # Max delta1_minority
  )

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

  # Generate initial values
  params_init <- foreach(minor = minority_groups,
                         .combine = dplyr::bind_rows) %:%
    foreach(sgn = c(-1, 1), .combine = bind_rows) %dopar% {
      tag_ <- paste(minor, sgn, sep = "_")

      if (!is.null(optim_params)) {
        init_ <- optim_params %>%
          dplyr::filter(tag == tag_) %>%
          dplyr::mutate(pars = map(optim, "par")) %>%
          dplyr::pull("pars") %>%
          `[[`(1)
      } else {
        init_ <- runif(length(params_lower),
                       min = params_lower,
                       max = params_upper)
      }

      dplyr::tibble(tag = tag_,  params = list(init_))
    }


  optim_res <- foreach(minor = minority_groups,
                       .combine = dplyr::bind_rows) %:%
    foreach(sgn = c(-1, 1), .combine = bind_rows) %dopar% {
      tag_ <- paste(minor, sgn, sep = "_")
      pars <- params_init %>%
        dplyr::filter(tag == tag_) %>%
        dplyr::pull("params") %>%
        `[[`(1)
      tibble(minor = minor,
             tag = tag_,
             optim = list(optim(pars,
                                .get_optim_fn(u, sgn = sgn, compare = c(base_group, minor), tag = tag_),
                                lower = params_lower, upper = params_upper)))
    }

  check_path(OPTIM_RDS)
  write_rds(optim_res, OPTIM_RDS)
  wfit2 <- function(f, d, ...) u$fit2(f, d, w = weights)

  d <- u$data

  # Filter and refactor data by grouping variable, as necessary
  if (!is.null(compare)) {
    # Check that compare is well defined
    check_levels <- sapply(compare, function(x) x %in% levels(d[[u$grouping]]))
    if (!all(check_levels)) {
      stop("Groups specified in compare do not exist:\n\t",
           compare[!check_levels],
           "\nAvailable groups are:\n\t",
           paste0(levels(d[[u$grouping]]), collapse = ", "))
    }

    target_group_ind <- d[[u$grouping]] %in% compare

    d <- d[target_group_ind, ]
    d[[u$grouping]] <- forcats::fct_drop(d[[u$grouping]])
    d[[u$grouping]] <- forcats::fct_relevel(d[[u$grouping]], compare)
  }

  if (is.factor(d[[u$grouping]])) {
    message(sprintf("Comparing %s=%s against %s={%s}",
                    u$grouping,
                    levels(d[[u$grouping]])[1],
                    u$grouping,
                    paste0(levels(d[[u$grouping]])[-1], collapse = ", ")))
  }

  # Initialize columns required for sensitize()
  d$treat <- d[[u$treatment]]

  if (is.null(ptreat)) {
    d$p_trt <- d$ptrt__
  } else if (length(ptreat) == nrow(d) | length(ptreat) == 1) {
    d$p_trt <- ptreat
  } else {
    stop("Bad specification of argument: ptreat")
  }

  if (is.null(resp_ctl)) {
    d$resp_ctl <- d$resp_ctl__
  } else if (length(resp_ctl) == nrow(d) | length(resp_ctl) == 1) {
    d$resp_ctl <- resp_ctl
  } else {
    stop("Bad specification of argument: resp_ctl")
  }

  if (is.null(resp_trt)) {
    d$resp_trt <- d$resp_trt__
  } else if (length(resp_trt) == nrow(d) | length(resp_trt) == 1) {
    d$resp_trt <- resp_trt
  } else {
    stop("Bad specification of argument: resp_trt")
  }

  # Appropriately expand parameters
  qs <- .expand_params(d[[u$grouping]], q)
  dps <- .expand_params(d[[u$grouping]], dp)
  d0s <- .expand_params(d[[u$grouping]], d0)
  d1s <- .expand_params(d[[u$grouping]], d1)

  sens_df <- rnr::sensitize(d, q = qs, dp = dps, d0 = d0s, d1 = d1s,
                            debug = TRUE)

  df_ <- dplyr::bind_rows(sens_df %>% dplyr::mutate(u = 0),
                          sens_df %>% dplyr::mutate(u = 1))

  df_[[u$treatment]] <- ifelse(df_$u == 0, df_$ptrt_u0__, df_$ptrt_u1__)

  if (u$risk_col == "resp_ctl") {
    beta__ <- df_$beta_ctl__
    delta__ <- df_$d0
  } else if (u$risk_col == "resp_trt") {
    beta__ <- df_$beta_trt__
    delta__ <- df_$d1
  } else {
    stop("Misspecified risk_col in undi object.\n\t",
         "Expected either resp_ctl or resp_trt\n\t",
         "Got: ", u$risk_col)
  }
  df_$risk__ <- beta__ + df_$u * delta__

  weights <- ifelse(df_$u == 0, 1 - df_$q, df_$q)

  df_$weights <- weights

  coefs <- .pull_coefs(df_,
                       u$treatment,
                       u$grouping,
                       c("risk__", u$controls),
                       fun = wfit2)

  ret <- coefs[grepl(u$grouping, coefs$term), ]

  if (debug) {
    ret <- list(df_, ret)
  }

  ret
}


#' Generate subroutine for computing sensitized race coefficients for single
#' minority group v. whites
#'
#' @param u undi object
#' @param sgn sign integer to multiply on scalar return value; used for
#'       controlling max/min optimization
#' @param compare vector of length 2, specifying the two groups to compare
#' @param verbose whether or not to print debug messages
#'       (0 = none, 1 = results only, 2 = everything)
#' @param tag string to tag output with (usefull for parallel output)
#'
#' @return
#'   Function that will return coefficient on minority group
.get_optim_fn <- function (u, sgn, compare, verbose = TRUE, tag = "fit") {
  function(params) {
    # Validate input
    if (length(compare) != 2) {
      stop("Can only get optim fn for comparison of two groups, got ",
           length(compare))
    }
    # Unpack parameters
    p <- .extract_params(params)
    qb  <- p$qb
    qm  <- p$qm
    ab  <- p$ab
    am  <- p$am
    d0b <- p$d0b
    d0m <- p$d0m
    d1b <- p$d1b
    d1m <- p$d1m

    if (verbose >= 2) {
      cat(sprintf(paste("%s: q=%.2f/%.2f, e(a)=%.2f/%.2f,",
                        "e(d0)=%.2f/%.2f, e(d1)=%.2f/%.2f\n"),
                  tag, qb, qm, exp(ab), exp(am), exp(d0b), exp(d0m),
                  exp(d1b), exp(d1m)))
    }

    ret <- undisens(u, compare = compare,
                    q = c(qb, qm),
                    dp = c(ab, am),
                    d0 = c(d0b, d0m),
                    d1 = c(d1b, d1m))
    ret <- ret[["estimate"]]

    ret <- ret * sgn

    if (verbose) {
      cat(sprintf("  %s coef: %.4f\n", tag, ret))
    }

    ret
  }
}

