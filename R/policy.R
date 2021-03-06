#' Create policy object
#'
#' @param formula a formula in the form of \code{treatment ~ grouping_variable +
#'   other predictors} where the LHS is the treatment column, first element on
#'   the RHS is the grouping variable (e.g., Race), and the remainder of the RHS
#'   specifies the predictors (features) for treatment/response models
#' @param data data frame to use; must include all the columns specified in
#'   \code{formula} and given in the \code{outcome} parameter
#' @param outcome name of outcome column in data
#' @param train either (1) a value between 0 and 1 representing the proportion
#'   of data to use in training, (2) the name of a column of characters "train"
#'   and "test" within \code{data} to use in splitting the data, or (3) a
#'   logical vector of equal length as \code{nrow(data)} used to index training
#'   data
#' @param risk One of \code{"resp_ctl"} or \code{"resp_trt"}, indicating which
#'   treatment regime should be used as the risk score (default:
#'   \code{resp_trt})
#' @param model character of modelling method to use for risk models; run
#'   \code{names(undi::models())} to get list of available models
#' @param down_sample (Optional) proportion (between 0 and 1) or number (greater
#'   than 1) of rows to sample, if down sampling the data; default is 1 (i.e.,
#'   use all data)
#' @param fit1 a function of the form f(formula, data, ...) used for fitting the
#'   first-stage model; overides \code{model}
#' @param pred1 a function of the form f(model, data, formula) used for
#'   generating predictions from the first-stage model; the formula argument can
#'   be ignored within the function body, but the function should still accept
#'   it; some prediction functions (e.g., glmnet) require the original formula;
#'   predictions should be on probability scale, while "risk" will always be on
#'   logit scale; overrides \code{model}
#' @param fit_ptreat a function of the form f(formula, data, ...) used for
#'   fitting propensity (probability of treatment) models. If not specified,
#'   \code{fit1} is used by default, with the provided \code{formula} argument.
#' @param pred_ptreat a function of the form f(model, data, formula) used for
#'   generating propensity predictions. If not specified, \code{pred1} is used
#'   by default.
#' @param ptreat (Optional) default value for probability of treatment;
#' @param resp_ctl (Optional)
#' @param resp_trt (Optional) default value for probability of response = 1
#'   given each treatment regime (\code{ctl}, \code{trt}); useful for cases
#'   where outcome under certain treatment regimes is deterministic (e.g.,
#'   probability of finding illegal weapon if NOT frisked is 0)
#' @param calibrate whether or not to use platt scaling to calibrate predictions
#' @param save_models whether or not fitted models should be returned
#' @param recast_df (Optional) whether \code{data} should be re-cast as a
#'   \code{data.frame}; this might be useful if the modeling function requires a
#'   model matrix, and \code{data} contains character columns --- recasting as a
#'   data frame will convert all character columns to factors
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
#'   \item{fit_ptreat}{function used to fit model for treatment propensity}
#'   \item{pred_ptreat}{function used to generate predictions for treatment
#'   propensity}\item{m_*}{if \code{save_models = TRUE}, each of the fitted
#'   models; otherwise set to \code{NULL}}\item{call}{the original call}
#'
#' @export
policy <-
  function(formula,
           data,
           outcome,
           train = 0.5,
           risk = "resp_trt",
           model = names(models()),
           down_sample = 1,
           fit1 = NULL,
           pred1 = NULL,
           fit_ptreat = NULL,
           pred_ptreat = NULL,
           ptreat = NULL,
           resp_ctl = NULL,
           resp_trt = NULL,
           calibrate = FALSE,
           save_models = FALSE,
           recast_df = FALSE,
           seed = round(stats::runif(1)*1e4),
           ...) {
    set.seed(seed)

    model <- match.arg(model)

    if (recast_df) {
      data <- as.data.frame(unclass(data))
    }

    data <- .down_sample(data, down_sample)

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
      fit1 <- models()[[model]]$fit
    }

    if (is.null(pred1)) {
      pred1 <- models()[[model]]$pred
    }

    if (is.null(fit_ptreat)) {
      fit_ptreat <- fit1
    }

    if (is.null(pred_ptreat)) {
      pred_ptreat <- pred1
    }

    # Split the data randomly, or by indexing vector, or use a predefined
    # column, based on the length of p_train
    if (length(train) == 1) {
      if(train > 0 & train < 1) {

        # Randomly assigned train/test splits stratified by treatment, group and outcome
        data <- data %>%
          group_by(!!as.name(treatment),
                   !!as.name(grouping),
                   !!as.name(outcome)) %>%
          mutate(fold__ = sample(c(
            rep('train', ceiling(train * n())),
            rep('test', n() - ceiling(train * n()))
          ),
          size = n())) %>%
          ungroup()

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
                m_resp_ctl = NULL,
                m_resp_trt = NULL,
                m_ptrt = NULL,
                ptreat = ptreat,
                resp_ctl = resp_ctl,
                resp_trt = resp_trt,
                calibrate = calibrate,
                save_models = save_models,
                call = match.call(expand.dots = F))
    class(ret) <- c("policy", class(ret))

    # Add policy estimates
    ret <- estimate_policy(ret, ...)

    return(ret)
  }

#' (Re)Estimate treatment and response surface models for a policy
#'
#' @param pol a policy object
#' @param features (Optional) character vector of features to use; if
#'   \code{NULL}, \code{pol$features} is used
#' @param ptreat (Optional) default value for probability of treatment;
#' @param resp_ctl (Optional)
#' @param resp_trt (Optional) default value for probability of response = 1
#'   given each treatment regime (\code{ctl}, \code{trt}); useful for cases
#'   where outcome under certain treatment regimes is deterministic (e.g.,
#'   probability of finding illegal weapon if NOT frisked is 0)
#' @param calibrate whether or not to use platt scaling to calibrate predictions
#' @param save_models whether or not fitted models should be returned
#' @param ... additional arguments  passed to model fitting functions \code{fit1}
#'   and \code{fit_ptreat}. See \code{\link{policy}}
#'
#' @return object of class \code{policy} where the relative columns in
#'   \code{$data} are updated with new predictions
#'
#' @export
estimate_policy <- function(pol,
                            features = NULL,
                            ptreat = pol$ptreat,
                            resp_ctl = pol$resp_ctl,
                            resp_trt = pol$resp_trt,
                            calibrate = pol$calibrate,
                            save_models = pol$save_models,
                            ...) {

    # Input validation
    if (!("policy" %in% class(pol))) {
      stop("Expected object of class policy")
    }

    # Get model fitting arguments from pol
    fit_args = overwrite_list(defaults = list(), pol$call$..., ...)

    d <- pol$data

    if (!is.null(features)) {
      pol$features <- features
    }

    train_ind <- d$fold__ == "train"
    ctl_train_ind <- train_ind & (d[[pol$treatment]] == 0)
    trt_train_ind <- train_ind & (d[[pol$treatment]] == 1)

    # Generate formulas for first and second stage models
    out_formula <- .make_formula(pol$outcome, c(pol$grouping, pol$features))
    trt_formula <- .make_formula(pol$treatment, c(pol$grouping, pol$features))

    # Fit first-stage models
    if (is.null(resp_trt)) {
      m1_trt <- do.call(pol$fit1, c(list(out_formula, d[trt_train_ind, ]), fit_args))

      d$resp_trt_pre_calib__ <- pol$pred1(m1_trt, d, out_formula)

      if (calibrate) {
        calib_formula <- .make_formula(pol$outcome, c("resp_trt_pre_calib__", pol$grouping))

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
      m1_ctl <- do.call(pol$fit1, c(list(out_formula, d[ctl_train_ind, ]), fit_args))

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
      m_ptrt <- do.call(pol$fit_ptreat, c(list(trt_formula, d[train_ind, ]), fit_args))

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

    # TODO: Check for possible model misspecification--
    # AUC < 0.6, miscalibration by group, few unique prediction
    # values, etc.

    pol$data <- d

    pol$data$risk__ <- .get_risk_col(pol)

    if (save_models) {
      pol$m_ptrt <- m_ptrt
      pol$m_resp_ctl <- m1_ctl
      pol$m_resp_trt <- m1_trt
    }

    return(pol)
}

