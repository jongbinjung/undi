#' Generate synthetic sample of a policy
#'
#' @param pol policy object
#' @param seed (Optional) random seed to set before generating synthetic
#'   observations
#'
#' @return new policy object with the treatment and outcome columns in
#'   \code{$data} updated with synthetic samples based on estimated values of
#'   \code{ptrt__}, \code{resp_ctl__}, and \code{resp_trt__}
#' @export
synthesize <- function(pol, seed = round(stats::runif(1)*1e4)) {
  ret <- pol

  N <- nrow(ret$data)

  set.seed(seed)

  ret$data[[ret$treatment]] <- stats::runif(N, 0, 1) <= ret$data$ptrt__

  resp_ctl <- stats::runif(N, 0, 1) <= ret$data$resp_ctl__
  resp_trt <- stats::runif(N, 0, 1) <= ret$data$resp_trt__

  ret$data[[ret$outcome]] <-
    ifelse(ret$data[[ret$treatment]], resp_trt, resp_ctl)

  return(ret)
}
