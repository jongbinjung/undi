#' Create policy object
#'
#' @param formula a formula in the form of \code{treatment ~ grouping_variable +
#'   other predictors} where the LHS is the treatment column, first element on
#'   the RHS is the grouping variable (e.g., Race), and the remainder of the RHS
#'   specifies the predictors for first-stage model
#' @param data data frame to use; must include all the columns specified in
#'   \code{formula} and given in the \code{outcome} parameter
#' @param outcome name of outcome column in data
#' @param train either (1) a value between 0 and 1 representing the proportion
#'   of data to use in training, (2) the name of a column of characters "train"
#'   and "test" within \code{data} to use in splitting the data, or (3) a
#'   logical vector of equal length as \code{nrow(data)} used to index training
#'   data
#' @param fit1 a function of the form f(formula, data, ...) used for fitting the
#'   first-stage model; using \code{gbm} by default
#' @param pred1 a function of the form f(model, data, formula) used for
#'   generating predictions from the first-stage model; the formula argument can
#'   be ignored within the function body, but the function should still accept
#'   it; some prediction functions (e.g., glmnet) require the original formula;
#'   predictions should be on probability scale, while "risk" will always be on
#'   logit scale
#' @param fit2 a function of the form f(formula, data, weights = NULL) used for
#'   fitting the second-stage model; using \code{glm} with \code{family =
#'   quasibinomial} by default; the \code{weights} argument is only used for
#'   sensitivity and *must* be initialized to \code{NULL} (or the equivalent of
#'   non-weighted fitting)
#' @param fit_ptreat a function of the form f(formula, data, ...) used for
#'   fitting propensity (probability of treatment) models. If not specified,
#'   \code{fit1} is used by default, with the provided \code{formula} argument.
#' @param pred_ptreat a function of the form f(model, data, formula) used for
#'   generating propensity predictions. If not specified, \code{pred1} is used
#'   by default.
#' @param risk One of \code{"resp_ctl"} or \code{"resp_trt"}, indicating which
#'   treatment regime should be used as the risk score (default:
#'   \code{resp_trt})
#' @param ptreat (Optional) default value for probability of treatment; if
#'   provided, it will override \code{fit_ptreat} and \code{pred_ptreat}
#' @param resp_ctl (Optional)
#' @param resp_trt (Optional) default value for probability of response = 1
#'   given each treatment regime (\code{ctl}, \code{trt}); useful for cases
#'   where outcome under certain treatment regimes is deterministic (e.g.,
#'   probability of finding illegal weapon if NOT frisked is 0)
#' @param calibrate whether or not to use platt scaling to calibrate predictions
#' @param save_models whether or not fitted models should be returned
#' @param seed random seed to use
#' @param ... additional arguments passed to first-stage model fitting function,
#'   \code{fit1} and \code{fit_ptreat}
#'
#' @return policy object \item{data}{original data frame, augmented with columns
#'   \code{fold__}, \code{ptrt__}, \code{resp_ctl__}, \code{resp_trt__}, and
#'   \code{risk__}, which contain train/test fold indicators for the first
#'   stage, treatment propesinty, first stage model probability predictions
#'   under control and treatment, and appropriate risk measure on logit-scale,
#'   respectively} \item{risk_col}{either "resp_ctl" or "resp_trt", indicating
#'   which was used for "risk"} \item{treatment}{name of column from \code{data}
#'   used as treatment indicator} \item{outcome}{name of column from \code{data}
#'   used as outcome indicator} \item{grouping}{name of column from \code{data}
#'   used as grouping variable} \item{features}{additional features used in
#'   first stage model} \item{fit1}{function used to fit first stage model}
#'   \item{pred1}{function used to generate predictions from first stage model}
#'   \item{fit2}{function used to fit second stage model}
#'   \item{fit_ptreat}{function used to fit model for treatment propensity}
#'   \item{pred_ptreat}{function used to generate predictions for treatment
#'   propensity}\item{m_*}{if \code{save_models = TRUE}, each of the fitted
#'   models; otherwise set to \code{NULL}}
#'
#' @export
policy <-
  function(formula,
           data,
           outcome,
           train = 0.5,
           fit1 = NULL,
           pred1 = NULL,
           fit2 = NULL,
           fit_ptreat = NULL,
           pred_ptreat = NULL,
           risk = "resp_trt",
           ptreat = NULL,
           resp_ctl = NULL,
           resp_trt = NULL,
           calibrate = FALSE,
           save_models = FALSE,
           seed = 1234,  # TODO: set seed randomly
           ...) {
    set.seed(seed)

    # Input validation
    if (length(risk) != 1 || !(risk %in% c("resp_ctl", "resp_trt"))) {
      stop("risk argument must be length 1 of either \"resp_ctl\" or \"resp_trt\"")
    }

    # Extract treatment/grouping/controls variables from formula
    features <- .extract_features(formula)

    treatment <- features$treat
    grouping <- features$group
    basefeats <- features$feats

    # TODO(jongbin): Sanity check columns in data

    # Check and initialize first/second-stage modeling functions
    if (is.null(fit1)) {
      fit1 <- function(f, d, ...) gbm::gbm(f, data = d, ...)
    }

    if (is.null(pred1)) {
      pred1 <- function(m, d, f)
        gbm::predict.gbm(m, d,
                         gbm::gbm.perf(m, plot.it = FALSE),
                         type = "response")
    }

    if (is.null(fit_ptreat)) {
      fit_ptreat <- fit1
    }

    if (is.null(pred_ptreat)) {
      pred_ptreat <- pred1
    }

    if (is.null(fit2)) {
      fit2 <-
        function(f, d, w = NULL) {
          if (is.null(w)) {
            stats::glm(f, d, family = stats::quasibinomial)
          } else {
            d$w <- w
            stats::glm(f, d, weights = w, family = stats::quasibinomial)
          }
        }
    }

    # Split the data randomly, or by indexing vector, or use a predefined
    # column, based on the length of p_train
    if (length(train) == 1) {
      if(train > 0 & train < 1) {
        # TODO: implement better random split
        data$fold__ <- sample(c("train", "test"),
                              size = nrow(data),
                              replace = TRUE,
                              prob = c(train, 1 - train))
      } else if (is.character(train)) {
        # TODO(jongbin): Check that specified column is "proper"
        data$fold__ <- data[[train]]
      } else {
        stop("train should either be between 0 and 1, or column name")
      }
    } else if (length(train) == nrow(data)) {
      data$fold__ <- ifelse(train, "train", "test")
    } else {
      stop("Wrong specification of argument train; see ?policy")
    }

    ret <- list(data = data,
                risk_col = risk,
                treatment = treatment,
                outcome = outcome,
                grouping = grouping,
                features = basefeats,
                fit1 = fit1,
                pred1 = pred1,
                fit_ptreat = fit_ptreat,
                pred_ptreat = pred_ptreat,
                fit2 = fit2,
                m_resp_ctl = NULL,
                m_resp_trt = NULL,
                m_ptrt = NULL)
    class(ret) <- c("policy", class(ret))

    # Add policy estimates
    ret <- estimate_policy(ret,
                           ptreat = ptreat,
                           resp_ctl = resp_ctl,
                           resp_trt = resp_trt,
                           calibrate = calibrate,
                           save_models = save_models)

    return(ret)
  }

#' (Re)Estimate treatment and response surface models for a policy
#'
#' @param pol a policy object
#' @param ptreat (Optional) default value for probability of treatment; if
#'   provided, it will override \code{fit_ptreat} and \code{pred_ptreat}
#' @param resp_ctl (Optional)
#' @param resp_trt (Optional) default value for probability of response = 1
#'   given each treatment regime (\code{ctl}, \code{trt}); useful for cases
#'   where outcome under certain treatment regimes is deterministic (e.g.,
#'   probability of finding illegal weapon if NOT frisked is 0)
#' @param calibrate whether or not to use platt scaling to calibrate predictions
#' @param save_models whether or not fitted models should be returned
#' @param ... additional arguments passed to model fitting functions \code{fit1}
#'   and \code{fit_ptreat}
#'
#' @return object of class \code{policy} where the relative columns in
#'   \code{$data} are updated with new predictions
#'
#' @export
estimate_policy <- function(pol,
                            ptreat = NULL,
                            resp_ctl = NULL,
                            resp_trt = NULL,
                            calibrate = FALSE,
                            save_models = FALSE,
                            ...) {
    # Input validation
    if (!("policy" %in% class(pol))) {
      stop("Expected object of class policy")
    }

    d <- pol$data

    train_ind <- d$fold__ == "train"
    ctl_train_ind <- train_ind & (d[[pol$treatment]] == 0)
    trt_train_ind <- train_ind & (d[[pol$treatment]] == 1)

    # Generate formulas for first and second stage models
    out_formula <- .make_formula(pol$outcome, c(pol$grouping, pol$features))
    trt_formula <- .make_formula(pol$treatment, c(pol$grouping, pol$features))

    # Fit first-stage models
    if (is.null(resp_trt)) {
      m1_trt <- pol$fit1(out_formula, d[trt_train_ind, ], ...)

      d$resp_trt_pre_calib__ <- pol$pred1(m1_trt, d, out_formula)

      if (calibrate) {
        calib_formula <- .make_formula(pol$outcome, "resp_trt_pre_calib__")

        d$resp_trt__ <-
          stats::predict.glm(stats::glm(calib_formula,
                                        data = d[trt_train_ind,],
                                        family = "binomial"),
                             d,
                             type = "response")
      } else {
        d$resp_trt__ <- d$resp_trt_pre_calib__
      }
    } else if (length(resp_trt) == nrow(d) | length(resp_trt) == 1) {
      message("Using custom values provided for resp_trt")
      m1_trt <- "Custom values of resp_trt provided"
      d$resp_trt__ <- resp_trt
    } else {
      stop("Bad specification of argument: resp_trt")
    }

    if (is.null(resp_ctl)) {
      m1_ctl <- pol$fit1(out_formula, d[ctl_train_ind, ], ...)

      d$resp_ctl_pre_calib__ <- pol$pred1(m1_ctl, d, out_formula)

      if (calibrate) {
        calib_formula <- .make_formula(pol$outcome, "resp_ctl_pre_calib__")

        d$resp_ctl__ <-
          stats::predict.glm(stats::glm(calib_formula,
                                        data = d[ctl_train_ind,],
                                        family = "binomial"),
                             d,
                             type = "response")
      } else {
        d$resp_ctl__ <- d$resp_ctl_pre_calib__
      }
    } else if (length(resp_ctl) == nrow(d) | length(resp_ctl) == 1) {
      message("Using custom values provided for resp_ctl")
      m1_ctl <- "Custom values of resp_ctl provided"
      d$resp_ctl__ <- resp_ctl
    } else {
      stop("Bad specification of argument: resp_ctl")
    }

    if (is.null(ptreat)) {
      m_ptrt <- pol$fit_ptreat(trt_formula, d[train_ind, ], ...)

      d$ptrt_pre_calib__ <- pol$pred_ptreat(m_ptrt, d, trt_formula)

      if (calibrate) {
        calib_formula <- .make_formula(pol$treatment, "ptrt_pre_calib__")

        d$ptrt__ <-
          stats::predict.glm(stats::glm(calib_formula,
                                        data = d[train_ind, ],
                                        family = "binomial"),
                             d,
                             type = "response")
      } else {
        d$ptrt__ <- d$ptrt_pre_calib__
      }
    } else if (length(ptreat) == nrow(d) | length(ptreat) == 1) {
      # TODO(jongbin): Provide warning for cases when fit/pred_ptreat is
      # specified but ignored
      message("Using custom values provided for ptreat")
      m_ptrt <- "Custom values of ptreat provided"
      d$ptrt__ <- ptreat
    } else {
      stop("Bad specification of argument: ptreat")
    }

    pol$data <- d

    pol$data$risk__ <- .get_risk_col(pol)

    if (save_models) {
      pol$m_ptrt <- m_ptrt
      pol$m_resp_ctl <- m1_ctl
      pol$m_resp_trt <- m1_trt
    }

    return(pol)
}
