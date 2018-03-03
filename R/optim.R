#' Find parameter values that min/maximize undi results
#'
#' Within specified range of sensitivity parameters, find the ones that achieve
#' minimum/maximum undi results
#'
#' @param r object of class \code{undi}
#' @param range_q 2D vector specifying min/max value of p(u = 1 | x)
#' @param range_dp 2D vector specifying min/max value of change in log-odds of
#'   treat = 1 if u = 1
#' @param range_d0 2D vector specifying min/max value of change in log-odds of
#'   response = 1 if treat = 0 and u = 1
#' @param range_d1 2D vector specifying min/max value of change in log-odds of
#'   response = 1 of treat = 1 and u = 1
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
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
  function(r,
           range_q = c(0, 1),
           range_dp = c(0, log(2)),
           range_d0 = c(0, log(2)),
           range_d1 = c(0, log(2)),
           base_group = NULL,
           minority_groups = NULL,
           debug = FALSE) {
  # Input validation
  if (!("undi" %in% class(r))) {
    stop("Expected object of class undi")
  }

  if (length(base_group) > 1) {
    stop("Specify a single base group.\n\tGot: ", base_group)
  }

  # range_q = c(0, 1)
  # range_dp = c(0, log(2))
  # range_d0 = c(0, log(2))
  # range_d1 = c(0, log(2))

  group_col <- r$data[[r$grouping]]
  members <- unique(group_col)

  check_groups <- sapply(c(base_group, minority_groups),
                         function(x) x %in% members)
  if (!all(check_groups)) {
    stop(sprintf("%s - not members of %s",
                 paste0(c(base_group, minority_groups)[!check_groups],
                        collapse = ","),
                 r$grouping))
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
  # TODO(jongbin): Allow the user to specify initial values?
  params_init <- foreach(minor = minority_groups,
                         .combine = dplyr::bind_rows) %:%
    foreach(sgn = c(-1, 1), .combine = dplyr::bind_rows) %dopar% {
      tag_ <- paste(minor, ifelse(sgn > 0, "min", "max"), sep = "_")

      init_ <- stats::runif(length(params_lower),
                            min = params_lower,
                            max = params_upper)

      dplyr::tibble(tag = tag_,  params = list(init_))
  }

  # Optimize for min/max over each minority group
  optim_res <- foreach(minor = minority_groups,
                       .combine = dplyr::bind_rows) %:%
    foreach(sgn = c(-1, 1), .combine = dplyr::bind_rows) %dopar% {
      tag_ <- paste(minor, ifelse(sgn > 0, "min", "max"), sep = "_")
      pars <- params_init %>%
        dplyr::filter(tag == tag_) %>%
        dplyr::pull("params") %>%
        `[[`(1)

      dplyr::tibble(minor = minor,
                    tag = tag_,
                    optim = list(
                      stats::optim(
                        pars,
                        .get_optim_fn(
                          r,
                          sgn = sgn,
                          compare = c(base_group, minor),
                          tag = tag_
                        ),
                        lower = params_lower,
                        upper = params_upper
                      )
                    ))
    }

  # Extract final results from
  coefs <-
    foreach(ip = 1:nrow(optim_res), .combine = dplyr::bind_rows) %dopar% {
    minor <- optim_res[ip, ][["minor"]]
    tag_ <- optim_res[ip, ][["tag"]]
    params <- optim_res[ip, ] %>%
      dplyr::mutate(pars = purrr::map(optim, "par")) %>%
      dplyr::pull("pars") %>%
      `[[`(1)

    fn <- .get_optim_fn(r, sgn = sgn, compare = c(base_group, minor),
                        return_scalar = FALSE)
    ret <- fn(params)
    ret$tag <- tag_
    ret
    }

  list(restuls = coefs, optim = optim_res)
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
#' @param return_scalar logical, whether to return a single scalar
#'       values (TRUE) or to return the full result from \code{undisens}
#'
#' @return
#'   Function that will return coefficient on minority group
.get_optim_fn <- function (u, sgn, compare, verbose = TRUE, tag = "fit",
                           return_scalar = TRUE) {
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
                    d1 = c(d1b, d1m),
                    verbose = FALSE)

    if (return_scalar) {
      ret <- ret[["estimate"]]

      if (verbose) {
        cat(sprintf("  %s coef: %.4f\n", tag, ret))
      }

      ret <- ret * sgn
    }

    ret
  }
}

