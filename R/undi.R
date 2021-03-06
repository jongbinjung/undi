#' undi: Testing for unjustified disparate impact
#'
#' The \code{undi} package implements a test for unjustified disparate impact in
#' policy decisions.
#'
#' @docType package
#' @name undi
#' @import rnr
#' @import dplyr
#' @import ggplot2
#' @importFrom rlang !!
NULL

if(getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      ":=",
      "ab",
      "am",
      "d0",
      "d0b",
      "d0m",
      "d1",
      "d1b",
      "d1m",
      "dp",
      "maxparam",
      "maxq",
      "minq",
      "minor",
      "qb",
      "qm",
      "risk_bin",
      "N",
      "mresp",
      "pout",
      "id_sens__",
      "fold__",
      "ptrt__",
      "risk__",
      "weights__",
      "resp_ctl__",
      "resp_trt__",
      "group",
      "sgn",
      "bound",
      "tag",
      "type",
      "ip",
      "optim",
      "term",
      "estimate",
      "controls",
      "lb",
      "ub",
      "cilb",
      "ciub",
      "base_cilb",
      "base_ciub",
      "odds_ratio",
      "std.error",
      "std.error.naive"
    )
  )
