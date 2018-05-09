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

#' \code{print} method of object class \code{sens} (parent of
#' \code{\link{optimsens}}, and \code{\link{gridsens}})
#'
#' @param x object of class sens
#' @param ... other stuff (for S3 consistency)
#'
#' @export
print.sens <- function(x, ...) {
  tags <- x$results %>%
    group_by(term) %>%
    arrange(estimate, .by_group = TRUE) %>%
    pull("tag")

  tag_length <- max(vapply(tags, nchar, FUN.VALUE = 0)) + 1
  ests <- format(x$results$estimate)

  for (tag in tags) {
    ind_res <- x$results$tag == tag

    est <- ests[ind_res]
    par <- x$results[ind_res, ]$pars[[1]]
    cat(format(tag, justify = "r", width = tag_length),
        " = ", est, " from ", .format_pars(par), "\n")
  }
}

#' Format optimsens parameters to a single line for printing
#'
#' @param pars optimsens parameters
#'
#' @return formatted string for printing
.format_pars <- function(pars) {
  sq <- sprintf("q=%.2f/%.2f", pars[["qb"]], pars[["qm"]])

  w <- getOption("digits") + 3

  if (pars[["ab"]] == pars[["am"]]) {
    sa <- paste0("a=", format(pars[["ab"]], width = w))
  } else {
    sa <- paste0("a=", format(pars[["ab"]], width = w), "/",
                 format(pars[["am"]], width = w))
  }

  if (pars[["d0b"]] == pars[["d0m"]]) {
    sd0 <- paste0("d0=", format(pars[["d0b"]], width = w))
  } else {
    sd0 <- paste0("d0=", format(pars[["d0b"]], width = w), "/", format(pars[["d0m"]], width = w))
  }

  if (pars[["d1b"]] == pars[["d1m"]]) {
    sd1 <- paste0("d1=", format(pars[["d1b"]], width = w))
  } else {
    sd1 <- paste0("d1=", format(pars[["d1b"]], width = w), "/", format(pars[["d1m"]], width = w))
  }

  return(paste(sq, sa, sd0, sd1, sep = " | "))
}
