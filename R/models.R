#' List of paired fit/pred methods by type (glm/glmnet/gbm)
#'
#' @return a list with model types (e.g., glm/gbm), each with appropriate
#'   \code{$fit} and \code{$pred} functions
#' @export
models <- function() {
  list(
    gbm = list(
      fit = function(f, d, ...)
        gbm::gbm(f, data = d, ...),
      pred = function(m, d, f)
        gbm::predict.gbm(m, d, gbm::gbm.perf(m, plot.it = FALSE),
                         type = "response")
    ),
    glm = list(
      fit = function(f, d, ...)
        stats::glm(f, data = d, family = stats::quasibinomial, ...),
      pred = function(m, d, f)
        stats::predict.glm(m, newdata = d, type = "response")
    ),
    lasso = list(
      fit = function(f, d, ...)
        fit_glmnet(f, d, family = "binomial", alpha = 1, ...),
      pred = function(m, d, f)
        predict4mm(m, d, f, type = "response", s = "lambda.min")
    ),
    ridge = list(
      fit = function(f, d, ...)
        fit_glmnet(f, d, family = "binomial", alpha = 0, ...),
      pred = function(m, d, f)
        predict4mm(m, d, f, type = "response", s = "lambda.min")
    ),
    sgd = list(
      fit = function(f, d, ...)
        fit_sgd(f, d, model = "glm",
                model.control = list(family = "binomial"), ...),
      pred = function(m, d, f) predict4mm(m, d, f, type = "response")
    )
  )
}

#' Wrapper for fitting \code{\link[sgd]{sgd}} models with standardized design
#' matrix
#'
#' @param f formula
#' @param d data object
#' @param ... additional arguments for \code{\link[sgd]{sgd}}
#'
#' @return an \code{sgd} object
#' @export
fit_sgd <- function(f, d, ...) {
  dm <- .get_mm(f, d)

  Y <- dm$Y
  x <- dm$x

  scaled <- .safe_scale(x)
  x[, scaled$ind] <- scaled$scaled
  m <- sgd::sgd(x, Y, ...)
  m$coef[scaled$ind] <-
    (m$coef[scaled$ind] + attr(scaled$scaled, "scaled:center")) *
    attr(scaled$scaled, "scaled:scale")

  return(m)
}


#' Wrapper for fitting \code{\link[glmnet]{glmnet}} models from a formula
#'
#' @param f formula
#' @param d data object
#' @param ... additional arguments for \code{\link[glmnet]{cv.glmnet}}
#'
#' @return an \code{\link[glmnet]{cv.glmnet}} object
#' @export
fit_glmnet <- function(f, d, ...) {
  dm <- .get_mm(f, d)

  Y <- dm$Y
  x <- dm$x

  m <- glmnet::cv.glmnet(x, Y, ...)

  return(m)
}

#' Wrapper to compute predictions from a model matrix
#'
#' @param m model object
#' @param d data frame to generate predictions on
#' @param f model formula; used to transform \code{d} to a \code{design matrix}
#' @param ... other arguments to pass to the predict function
#'
#' @return generated predictions
#' @export
predict4mm <- function(m, d, f, ...) {
  mm <- stats::model.matrix(f, d)

  stats::predict(m, mm, ...)
}


#' Get model matrix and response given a formula and data frame
#'
#' @param f model formula
#' @param d data frame
#'
#' @return a list with \code{Y} and \code{x}, the response and design matrix,
#'   respectively
.get_mm <- function(f, d) {
  mf <- stats::model.frame(f, d)

  # Extract Y values
  Y <- stats::model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) {
      names(Y) <- nm
    }
  }

  # Build and standardize model matrix x
  mt <- attr(mf, "terms")
  if (!stats::is.empty.model(mt)) {
    x <- stats::model.matrix(mt, mf)
  } else {
    stop("No valid terms found in data for specified formula: ", format(f))
  }

  list(Y = Y, x = x)
}
