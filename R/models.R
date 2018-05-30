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
        undi::fit_glmnet(f, d, family = "binomial", alpha = 1, ...),
      pred = function(m, d, f)
        undi::predict4mm(m, d, f, type = "response", s = "lambda.min")
    ),
    ridge = list(
      fit = function(f, d, ...)
        undi::fit_glmnet(f, d, family = "binomial", alpha = 0, ...),
      pred = function(m, d, f)
        undi::predict4mm(m, d, f, type = "response", s = "lambda.min")
    ),
    sgd = list(
      fit = function(f, d, model.control = list(family = "binomial"),
                     sgd.control = list(lr = 'adagrad', reltol = 1e-8, shuffle = T), ...)
        undi::fit_sgd(f, d, model = "glm",
                      model.control = model.control,
                      sgd.control = sgd.control, ...),
      pred = function(m, d, f)
        undi::predict4mm(m, d, f, type = "response")
    )
  )
}


#' Given a policy and (optional) controls, generate a rad_control object
#'
#' @param pol a \code{\link{policy}} object
#' @param fit_fn string indicating the rad estimation model/procedure used.
#'   \code{*_coef} methods use models without interaction between risk and
#'   group, and return the coeficient on group membership. \code{*_avg} methods
#'   will fit more flexible models (possibly with interactions), and compute
#'   average ratios across the population. (TODO: better documentation is
#'   expected)
#' @param controls character vector of additional controls to consider in the
#'   second-stage model
#' @param use_speedglm whether or not to use \code{speedglm}, instead of
#'   \code{stats::glm}, in cases where N > 2P (see details)
#'
#' @details speedglm can potentially speed-up computation significantly, but
#'   only in cases where the number of rows is somewhat greater than the number
#'   of features (specifically, when N > 2P). In terms of FLOPs at each Fisher
#'   iteration, stats::glm requires (2np^2 - (2/3)p^3) FLOPS vs, (np^2 +
#'   (4/3)p^3) for speedglm.
#'
#' @return a \code{rad_control} object constructed of \item{formula}{the formula
#'   used in model fitting} \item{label}{a character label associated with the
#'   model fit type} \item{grouping}{column name of group, as specified in
#'   \code{pol$grouping}} \item{fit}{a function of the form f(d, w = NULL, ...)
#'   for fitting a model with training data \code{d}} \item{pred}{function of
#'   the form g(m, d) for generating predictions for data \code{d} with model
#'   \code{m}} \item{method}{character string describing the method to use}
#' @export
rad_control <-
  function(pol,
           fit_fn = c("logit_coef", "gam_coef", "decbin_coef",
                      "logit_avg", "gam_avg"),
           controls = NULL,
           use_speedglm = TRUE) {
  fit_fn <- match.arg(fit_fn)

  switch (fit_fn,
    logit_coef = {
      feats <- c(pol$grouping, "risk__", controls)

      f <- .make_formula(pol$treatment, feats)

      label <- paste(feats, collapse = ", ")

      ret <- list(
        formula = f,
        label = label,
        fit = NA,
        pred = function() stop("Never call pred() on a coef method for rad!"),
        method = "coef"
      )

      # FLOPs at each iteration:
      #   stats::glm: 2np^2 - (2/3)p^3
      #   speedglm: np^2 + (4/3)p^3
      # Use stats::glm only if n < 2p
      if (use_speedglm && nrow(pol$data) > 2 * length(feats)) {
        ret$fit <- function(d, w = NULL, ...) {
          if (is.null(w)) {
            speedglm::speedglm(f, d, family = stats::quasibinomial(), ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            speedglm::speedglm(f, d, weights = w,
                               family = stats::quasibinomial(),
                               ...)
          }
        }
      } else {
        ret$fit <- function(d, w = NULL, ...) {
          if (is.null(w)) {
            stats::glm(f, d, family = stats::quasibinomial(), ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            stats::glm(f, d, weights = w, family = stats::quasibinomial(), ...)
          }
        }
      }
    },
    gam_coef = {
      feats <- c(pol$grouping, "s(risk__)", controls)

      f <- .make_formula(pol$treatment, feats)

      label <- paste(feats, collapse = ", ")

      ret <- list(
        formula = f,
        label = label,
        fit = function(d, w = NULL, ...) {
          if (is.null(w)) {
            mgcv::gam(f, data = d, family = stats::quasibinomial(), ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            mgcv::gam(f, data = d, weights = w, family = stats::quasibinomial(), ...)
          }
        },
        pred = function() stop("Never call pred() on a coef method for rad!"),
        method = "coef"
      )
    },
    decbin_coef = {
      deciles <- pol$data %>%
        filter(fold__ == "test") %>%
        pull(risk__) %>%
        stats::quantile(seq(.1, .9, .1))


      feats <- c(pol$grouping, "riskbin__", controls)

      f <- .make_formula(pol$treatment, feats)

      label <- paste(feats, collapse = ", ")

      ret <- list(
        formula = f,
        label = label,
        fit = function(d, w = NULL, ...) {
          d$riskbin__ <- cut(d$risk__, c(-Inf, deciles, Inf))
          if (is.null(w)) {
            mgcv::gam(f, data = d, family = stats::quasibinomial(), ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            mgcv::gam(f, data = d, weights = w, family = stats::quasibinomial(), ...)
          }
        },
        pred = function() stop("Never call pred() on a coef method for rad!"),
        method = "coef"
      )

      if (use_speedglm && nrow(pol$data) > 2 * length(feats)) {
        ret$fit <- function(d, w = NULL, ...) {
          d$riskbin__ <- cut(d$risk__, c(-Inf, deciles, Inf))
          if (is.null(w)) {
            speedglm::speedglm(f, d, family = stats::quasibinomial(), ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            speedglm::speedglm(f, d, weights = w,
                               family = stats::quasibinomial(),
                               ...)
          }
        }
      } else {
        ret$fit <- function(d, w = NULL, ...) {
          d$riskbin__ <- cut(d$risk__, c(-Inf, deciles, Inf))
          if (is.null(w)) {
            stats::glm(f, d, family = stats::quasibinomial(), ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            stats::glm(f, d, weights = w, family = stats::quasibinomial(), ...)
          }
        }
      }
    },
    logit_avg = {
      feats <- c(pol$grouping, "risk__", controls)

      f <- .make_formula(pol$treatment,
                         c(feats, paste0(pol$grouping, ":", "risk__")))

      label <- paste(feats, collapse = ", ")

      ret <- list(
        formula = f,
        label = label,
        fit = function(d, w = NULL, ...) {
          if (is.null(w)) {
            stats::glm(f, d, family = stats::quasibinomial, ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            stats::glm(f, d, weights = w, family = stats::quasibinomial, ...)
          }
        },
        pred = function(m, d) {
          stats::predict.glm(object = m,
                             newdata = d,
                             type = "response")
        },
        method = "avg"
      )
    },
    gam_avg = {
      feats <- c(pol$grouping, "gam::s(risk__)", controls)

      f <- .make_formula(pol$treatment,
                         c(feats, paste0(pol$grouping, ":", "risk__")))

      label <- paste(feats, collapse = ", ")

      ret <- list(
        formula = f,
        label = label,
        fit = function(d, w = NULL, ...) {
          if (is.null(w)) {
            gam::gam(f, data = d, family = stats::quasibinomial, ...)
          } else {
            # Make sure that the weights exist in d at the time of call
            d$w <- w
            gam::gam(f, data = d, weights = w, family = stats::quasibinomial, ...)
          }
        },
        pred = function(m, d) {
          gam::predict.Gam(m, d, type = "response")
        },
        method = "avg"
      )
    },
    stop("Failed to construct disparate impact model functions\n",
         "\t(Unknown fit_fn?)")
  )

  ret$grouping <- pol$grouping

  class(ret) <- c("rad_control", class(ret))
  return(ret)
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

  x <- .safe_scale(x)
  m <- sgd::sgd(x, Y, ...)

  m$scale = attr(x, "scaled:scale")
  m$center = attr(x, "scaled:center")
  m$formula = f

  m$coefficients[1] <- m$coefficients[1] -
    sum(m$center[-1]*m$coefficients[-1]/m$scale[-1])

  m$coefficients[-1] <-
    m$coefficients[-1]/m$scale[-1]

  m$coefficients <- stats::setNames(m$coefficients[,1], colnames(x))

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

  preds <- stats::predict(m, mm, ...)

  if (is.null(dim(preds))) {
    return(preds)
  } else if (length(dim(preds)) == 1) {
    dim(preds) <- NULL
    return(preds)
  } else if (length(dim(preds)) == 2 & dim(preds)[2] == 1) {
    dim(preds) <- NULL
    return(preds)
  } else {
    stop("Predictions have unexpected dimensions for model class: ", class(m))
  }
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
