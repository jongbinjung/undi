#' \code{print} method of object class \code{\link{policy}}
#'
#' @param x an object of class policy
#' @param ... other stuff
#'
#' @export
print.policy <- function(x, ...) {
  # TODO: Implement me
  warning("TODO: print.policy() needs to be implemented properly")
  d <- x$data

  train_ind <- d$fold__ == "train"
  ctl_train <- train_ind & (d[[x$treatment]] == 0)
  trt_train <- train_ind & (d[[x$treatment]] == 1)

  ctl_test <- (!train_ind) & (d[[x$treatment]] == 0)
  trt_test <- (!train_ind) & (d[[x$treatment]] == 1)

  s <- paste("Models have train/test AUC",
             "\tptreat  : %s / %s",
             "\tresp_ctl: %s / %s",
             "\tresp_trt: %s / %s",
             "",
             sep = "\n")

  cat(
    sprintf(
      s,
      .compute_auc(d[train_ind, "ptrt__"], d[train_ind, ][[x$treatment]]),
      .compute_auc(d[!train_ind, "ptrt__"], d[!train_ind, ][[x$treatment]]),
      .compute_auc(d[ctl_train, "resp_ctl__"], d[ctl_train, ][[x$outcome]]),
      .compute_auc(d[ctl_test, "resp_ctl__"], d[ctl_test, ][[x$outcome]]),
      .compute_auc(d[trt_train, "resp_trt__"], d[trt_train, ][[x$outcome]]),
      .compute_auc(d[trt_test, "resp_trt__"], d[trt_test, ][[x$outcome]])
    )
  )
}
