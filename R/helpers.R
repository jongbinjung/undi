.make_formula <- function(y, vars) {
  stats::as.formula(paste(y, "~", paste(vars, collapse = "+")))
}

logit <- stats::binomial()$linkfun
inv_logit <- stats::binomial()$linkinv

#' Fit second stage model and extract coefficients
#'
#' Get tidy coefficients for cn_lhs ~ cn_tgt + controls given a data frame and
#' additional controls
#'
#' @param d data frame that has all columns referenced in cn_lhs, cn_tgt, and
#'   controls
#' @param cn_lhs name of column to use as LHS
#' @param cn_tgt the target column of interest
#' @param controls character vector of additional columns to control for
#' @param fun function (formula, data, ...) used to fit model; must return a
#'   model object that can be handled by broom::tidy() (default: glm)
#'
#' @return tidy dataframe of the glm model
.pull_coefs <-
  function(d,
           cn_lhs,
           cn_tgt,
           controls = NULL,
           fun = NULL) {
  if (is.null(fun)) {
    stop("Second-stage fitting function not specified")
  }

  f <- .make_formula(cn_lhs, c(controls, cn_tgt))
  lbl_controls <- ifelse(is.null(controls), "None",
                         ifelse(length(controls) > 5, "Kitchen sink",
                                paste(controls, collapse = ", ")))
  fun(f, d) %>%
    broom::tidy() %>%
    dplyr::mutate(controls = lbl_controls)
}

#' Extract LHS, first element of RHS, and remainder of RHS from a formula
#'
#' Used for extracting relevant column names for features/grouping
#'
#' @param formula a formula where the first element of the RHS represents a
#'   grouping variable
#'
#' @return list of elements treat, group, and feats
.extract_features <- function(formula) {
  treat <- as.character(lazyeval::f_lhs(formula))
  group <- labels(stats::terms(formula))[1]
  feats <- labels(stats::terms(formula))[-1]

  list(treat = treat, group = group, feats = feats)
}

#' Expand parameters by length
#'
#' Used for determining how to apply sensitivity parameters
#'
#' @param g vector of grouping variable
#' @param p paramter in size of either 1, equal to the number of levels in
#'   \code{g}, or equal to the length of \code{g}
#'
#' @return vector of parameter \code{p} appropriately expanded
.expand_params <- function(g, p) {
  len <- length(p)

  if (len == 1) {
    return(rep(p, length(g)))
  } else if (is.factor(g) & len == length(levels(g))) {
    return(p[as.numeric(g)])
  } else if (len == length(g)) {
    return(p)
  } else {
    stop("Mis-specified parameter assignments to grouping variable.\n\t",
         "Expecting parameters of length 1",
         if(is.factor(g)) paste0(", ", length(levels(g)), ","),
         " or ", length(g), "\n\t",
         "Received parameter of length: ", length(p))
  }
}

