#' Sensitivity to unobserved confounders with specified parameters
#'
#' @param pol object of class \code{policy}
#' @param q p(u = 1 | x) (see Details)
#' @param dp change in log-odds of treat = 1 if u = 1 (see Details)
#' @param d0 change in log-odds of response = 1 if treat = 0 and u = 1 (see
#'   Details)
#' @param d1 change in log-odds of response = 1 of treat = 1 and u = 1 (see
#'   Details)
#' @param compare (Optional) character vector of groups to compare; the data
#'   will be filtered such that only specified groups are compared, and the
#'   \code{grouping} variable will be refactored such that the levels preserve
#'   the specified order, e.g., \code{compare = c("white", "black")} will make
#'   \code{"white"} the base group
#' @param ptreat (Optional) default value for probability of treatment; if
#'   provided, it will override fitted values in \code{pol$data}
#' @param resp_ctl (Optional)
#' @param resp_trt (Optional) default value for probability of response = 1
#'   given each treatment regime (\code{ctl}, \code{trt}); useful for cases
#'   where outcome under certain treatment regimes is deterministic (e.g.,
#'   probability of finding illegal weapon if NOT frisked is 0); if provided, it
#'   will override fitted values in \code{pol$data}
#' @param controls vector of legitimate controls to use; the ones specified
#'   within the policy object will be used if not specified
#' @param naive_se logical flag, if TRUE, return std.error from naive
#'   (non-sensitivity) test, as well as std.errors from final weighted
#'   regression
#' @param fit_fn string indicating the fitting proceedure used.
#' @param verbose logical flag, if TRUE, print relevant messages for user
#' @param debug logical flag, if TRUE, returns a list of results and the
#'   expanded data frame used to fit model
#' @param ... additional arguments to pass to \code{fit} function from
#'   \code{\link{di_model}} for fine-tuning
#'
#' @details All sensitivity parameters (\code{q, dp, d0, d1}) can be provided in
#'   one of three formats, determined by the \code{length} of each argument:
#'   \describe{ \item{if \code{length(arg) = 1}}{single value applied to all
#'   observations (rows)} \item{if \code{length(arg) = }number of levels in
#'   grouping variable}{each parameter setting applied to corresponding level in
#'   group} \item{if \code{length(arg) = nrow(pol$data)}}{each parameter applied
#'   to corresponding rows}} Note that if \code{compare} is specified, the
#'   number of grouping levels is effectively the length of \code{compare}
#'
#' @return \code{tidy} dataframe of second-stage model coefficients after
#'   applying sensitivity parameters, with a nested column of sensitivity params
#'
#' @export
sensitivity <-
  function(pol,
           q,
           dp,
           d0,
           d1,
           compare = NULL,
           ptreat = NULL,
           resp_ctl = NULL,
           resp_trt = NULL,
           controls = NULL,
           naive_se = TRUE,
           fit_fn = c("logit"),
           verbose = interactive(),
           debug = FALSE,
           ...) {

  if (!("policy" %in% class(pol))) {
    stop("Expected object of class policy")
  }

  fit_fn <- match.arg(fit_fn)
  dm <- di_model(pol, controls, fit_fn = fit_fn, ...)

  if (is.null(controls)) {
    controls <- pol$controls
  }

  if (naive_se) {
    # Get standard errors from non-sensitivity adjusted regression
    od <- pol$data
    od <- od[od$fold__ == "test", ]

    if (!is.null(compare)) {
      # Check that compare is well defined
      check_levels <- sapply(compare,
                             function(x) x %in% levels(od[[pol$grouping]]))

      if (!all(check_levels)) {
        stop("Groups specified in compare do not exist:\n\t",
             compare[!check_levels],
             "\nAvailable groups are:\n\t",
             paste0(levels(od[[pol$grouping]]), collapse = ", "))
      }

      target_group_ind <- od[[pol$grouping]] %in% compare

      od <- od[target_group_ind, ]
      od[[pol$grouping]] <- forcats::fct_drop(od[[pol$grouping]])
      od[[pol$grouping]] <- forcats::fct_relevel(od[[pol$grouping]], compare)
    }

    naive_coefs <- .compute_estimate(
      d = od,
      cn_group = pol$grouping,
      dm
    )

    naive_coefs <- naive_coefs[, c("term", "std.error")] %>%
      dplyr::mutate(controls = paste(c(pol$grouping, "risk__", controls),
                                     collapse = ", "))
  }

  sens_pol <- sensitize(pol, q = q, dp = dp, d0 = d0, d1 = d1,
                        compare = compare,
                        ptreat = ptreat,
                        resp_ctl = resp_ctl,
                        resp_trt = resp_trt)

  d <- sens_pol$sens_data
  d <- d[d$fold__ == "test", ]

  if (is.factor(d[[pol$grouping]])) {
    if (verbose) {
      message(sprintf("Comparing %s=%s against %s={%s}",
                      pol$grouping,
                      levels(d[[pol$grouping]])[1],
                      pol$grouping,
                      paste0(levels(d[[pol$grouping]])[-1], collapse = ", ")))
    }
  }

  d_ <- dplyr::filter(d, weights__ != 0)

  coefs <- .compute_estimate(d = d_,
                             cn_group = pol$grouping,
                             dm = dm,
                             weighted = TRUE)

  # TODO(jongbin): "sgd" should be re-implemented or removed?
  # } else if (fit_fn == 'sgd') {
  #   form <- .make_formula(NULL, c('risk__', pol$grouping, controls))
  #   X <- stats::model.matrix(form, d_)
  #
  #   sgd_result <- fit_sgd(
  #     .make_formula(pol$treatment,
  #                   c('risk__', pol$grouping, controls)),
  #     d_,
  #     model = 'glm',
  #     model.control = list(family = 'binomial', weights = d_$weights),
  #     sgd.control = list(
  #       lr = 'adagrad',
  #       reltol = 1e-8,
  #       shuffle = TRUE
  #     )
  #   )
  #
  #   coefs <- data.frame(
  #     term = colnames(X),
  #     estimate = sgd_result$coefficients,
  #     std.error = NA,
  #     controls = paste(c(pol$grouping, 'risk__', controls), collapse = ", "))
  # }


  # Replace std.error with naive estimates
  coefs <- coefs[, c("term", "estimate", "std.error", "controls")]

  if (naive_se) {
    coefs <- merge(coefs,
                   naive_coefs,
                   by = c("term", "controls"),
                   suffixes = c(".weighted", ".naive"))
  }

  ret <- coefs[grepl(pol$grouping, coefs$term), ]

  if (debug) {
    ret <- list(d_, ret)
  }

  ret
  }


#' Compute the sensitivity-adjusted estimates of predicted outcome given
#' treatment/control for a policy object
#'
#' @param obj data to sensitize
#' @param q p(u = 1 | x) (see Details)
#' @param dp change in log-odds of treat = 1 if u = 1 (see Details)
#' @param d0 change in log-odds of response = 1 if treat = 0 and u = 1 (see
#'   Details)
#' @param d1 change in log-odds of response = 1 of treat = 1 and u = 1 (see
#'   Details)
#' @param compare (Optional) character vector of groups to compare; the data
#'   will be filtered such that only specified groups are compared, and the
#'   \code{grouping} variable will be refactored such that the levels preserve
#'   the specified order, e.g., \code{compare = c("white", "black")} will make
#'   \code{"white"} the base group
#' @param ptreat (Optional) default value for probability of treatment; if
#'   provided, it will override fitted values in \code{pol$data}
#' @param resp_ctl (Optional)
#' @param resp_trt (Optional) default value for probability of response = 1
#'   given each treatment regime (\code{ctl}, \code{trt}); useful for cases
#'   where outcome under certain treatment regimes is deterministic (e.g.,
#'   probability of finding illegal weapon if NOT frisked is 0); if provided, it
#'   will override fitted values in \code{pol$data}
#' @param ... additional arguments are ignored
#'
#' @details All sensitivity parameters (\code{q, dp, d0, d1}) can be provided in
#'   one of three formats, determined by the \code{length} of each argument:
#'   \describe{ \item{if \code{length(arg) = 1}}{single value applied to all
#'   observations (rows)} \item{if \code{length(arg) = }number of levels in
#'   grouping variable}{each parameter setting applied to corresponding level in
#'   group} \item{if \code{length(arg) = nrow(pol$data)}}{each parameter applied
#'   to corresponding rows}} Note that if \code{compare} is specified, the
#'   number of grouping levels is effectively the length of \code{compare}
#'
#' @return a new \code{sensitive_policy} object, which inherits \code{policy},
#'   with the data element updated according to sensitivity parameters; all
#'   other aspects of the original \code{policy} object are preserved.
#' @export
sensitize.policy <-
  function(obj,
           q,
           dp,
           d0,
           d1,
           compare = NULL,
           ptreat = NULL,
           resp_ctl = NULL,
           resp_trt = NULL,
           ...) {
  if (!("policy" %in% class(obj))) {
    stop("Expected object of class policy")
  }

  # Pull out data, and tag with internal ID for future consistency (e.g., plots)
  d <- obj$data
  d$id_sens__ <- 1:nrow(d)

  # Filter and refactor data by grouping variable, as necessary
  if (!is.null(compare)) {
    # Check that compare is well defined
    check_levels <- sapply(compare, function(x) x %in% levels(d[[obj$grouping]]))

    if (!all(check_levels)) {
      stop("Groups specified in compare do not exist:\n\t",
           compare[!check_levels],
           "\nAvailable groups are:\n\t",
           paste0(levels(d[[obj$grouping]]), collapse = ", "))
    }

    target_group_ind <- d[[obj$grouping]] %in% compare

    d <- d[target_group_ind, ]
    d[[obj$grouping]] <- forcats::fct_drop(d[[obj$grouping]])
    d[[obj$grouping]] <- forcats::fct_relevel(d[[obj$grouping]], compare)
  }

  # Initialize columns required for sensitize()
  d$treat <- d[[obj$treatment]]

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
  qs <- .expand_params(d[[obj$grouping]], q)
  dps <- .expand_params(d[[obj$grouping]], dp)
  d0s <- .expand_params(d[[obj$grouping]], d0)
  d1s <- .expand_params(d[[obj$grouping]], d1)

  new_d <- sensitize(d, q = qs, dp = dps, d0 = d0s, d1 = d1s, debug = TRUE)

  df_ <- dplyr::bind_rows(new_d %>% dplyr::mutate(u = 0),
                          new_d %>% dplyr::mutate(u = 1))

  df_[[obj$treatment]] <- ifelse(df_$u == 0, df_$ptrt_u0__, df_$ptrt_u1__)

  if (obj$risk_col == "resp_ctl") {
    beta__ <- df_$beta_ctl__
    delta__ <- df_$d0
  } else if (obj$risk_col == "resp_trt") {
    beta__ <- df_$beta_trt__
    delta__ <- df_$d1
  } else {
    stop("Misspecified risk_col in policy object.\n\t",
         "Expected either resp_ctl or resp_trt\n\t",
         "Got: ", obj$risk_col)
  }

  df_$risk__ <- beta__ + df_$u * delta__

  df_$weights__ <- ifelse(df_$u == 0, 1 - df_$q, df_$q)

  # Update and reclass policy object
  obj$data <- d
  obj$sens_data <- df_
  obj$params <- list(q = q, d = dp, d0 = d0, d1 = d1, compare = compare)

  class(obj) <- c("sensitive_policy", class(obj))

  return(obj)
}
