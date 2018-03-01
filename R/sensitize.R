#' Sensitivity to unobserved confounders with specified parameters
#'
#' @param u object of class \code{undi}
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
#'   provided, it will override fitted values in \code{u$data}
#' @param resp_ctl (Optional)
#' @param resp_trt (Optional) default value for probability of response = 1
#'   given each treatment regime (\code{ctl}, \code{trt}); useful for cases
#'   where outcome under certain treatment regimes is deterministic (e.g.,
#'   probability of finding illegal weapon if NOT frisked is 0); if provided, it
#'   will override fitted values in \code{u$data}
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
undisens <-
  function(u,
           q,
           dp,
           d0,
           d1,
           compare = NULL,
           ptreat = NULL,
           resp_ctl = NULL,
           resp_trt = NULL,
           debug = FALSE) {

  if (!("undi" %in% class(u))) {
    stop("Expected object of class undi")
  }

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
}
