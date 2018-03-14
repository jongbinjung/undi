#' Summary of optimsens object
#'
#' @param x optimsens object
#' @param ... other stuff
#'
summary.optimsens <- function(x, ...) {
  format(x$results %>%
    tidyr::separate(tag, into = c("group", "bound"), sep = "_") %>%
    dplyr::select(term, estimate, bound) %>%
    tidyr::spread(bound, estimate) %>%
    dplyr::select(term, min, max))
}
