#' undi: Testing for unjustified disparate impact
#'
#' The \code{undi} package implements a test for unjustified disparate impact in
#' policy decisions.
#'
#' @docType package
#' @name undi
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom foreach %:%
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
NULL

if(getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      "sgn",
      "tag",
      "ip",
      "optim",
      "term",
      "estimate",
      "lb",
      "ub",
      "cilb",
      "ciub",
      "odds_ratio",
      "std.error.naive"
    )

  )
