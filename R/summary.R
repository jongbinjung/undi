#' Summary of sens object
#'
#' @param object \code{sens} object
#' @param ... other stuff
#'
#' @return summary data frame with min/max values of each group term
#'
#' @export
summary.sens <- function(object, ...) {
  object$results %>%
    tidyr::separate(tag, into = c("group", "bound"), sep = "_") %>%
    dplyr::select(term, estimate, bound) %>%
    tidyr::spread(bound, estimate) %>%
    dplyr::select(term, min, max)
}
