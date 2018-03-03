#' undi: Testing for unjustified disparate impact
#'
#' The \code{undi} package implements a test for unjustified disparate impact in
#' policy decisions.
#'
#' @docType package
#' @name undi
#' @importFrom dplyr %>%
#' @importFrom foreach %:%
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
NULL

if(getRversion() >= "2.15.1")
  utils::globalVariables(c("sgn", "tag", "ip", "optim")
  )
