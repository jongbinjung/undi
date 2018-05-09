# Implementation of common (dplyr) verbs for policy manipulation

#' Generate a new policy after applying various dplyr verbs
#'
#' @param .data a policy object
#' @param x a policy object
#' @param ... arguments passed to corresponding dplyr verb
#'
#' @details Note that the policy object is \emph{never} refit, and only the
#'   \code{$data} element is manipulated. In certain cases (e.g., after a
#'   filter), the policy might need to be refit via
#'   \code{\link{estimate_policy}}.
#'
#' @return new policy object with altered data
#' @name dplyr.policy
NULL

#' @rdname dplyr.policy
#' @export
filter.policy <- function(.data, ...) {
  .data$data <- filter(.data$data, ...)

  return(.data)
}

#' @rdname dplyr.policy
#' @export
mutate.policy <- function(.data, ...) {
  .data$data <- mutate(.data$data, ...)

  return(.data)
}

#' @rdname dplyr.policy
#' @export
group_by.policy <- function(.data, ...) {
  .data$data <- group_by(.data$data, ...)

  return(.data)
}

#' @rdname dplyr.policy
#' @export
ungroup.policy <- function(x, ...) {
  x$data <- ungroup(x$data, ...)

  return(x)
}

