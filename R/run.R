#' Create policy object
#'
#' @param formula a formula in the form of \code{treatment ~ grouping_variable +
#'   other predictors} where the LHS is the treatment column, first element on
#'   the RHS is the grouping variable (e.g., Race), and the remainder of the RHS
#'   specifies the predictors for first-stage model
#' @param data data frame to use; must include all the columns specified in
#'   \code{formula} and given in the \code{outcome} parameter
#' @param outcome name of outcome column in data
#' @param controls character vector of additional controls to consider in the
#'   second-stage model
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
#'   which was used for "risk"}\item{m1_*}{fitted first stage model for response
#'   give ctl/trt} \item{treatment}{name of column from \code{data} used as
#'   treatment indicator} \item{outcome}{name of column from \code{data} used as
#'   outcome indicator} \item{grouping}{name of column from \code{data} used as
#'   grouping variable} \item{features}{additional features used in first stage
#'   model} \item{controls}{legitimate controls used in second stage model}
#'   \item{fit1}{function used to fit first stage model} \item{pred1}{function
#'   used to generate predictions from first stage model} \item{fit2}{function
#'   used to fit second stage model} \item{fit_ptreat}{function used to fit
#'   model for treatment propensity} \item{pred_ptreat}{function used to
#'   generate predictions for treatment propensity}
#'
#' @export
policy <-
  function(formula,
           data,
           outcome,
           controls = NULL,
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

    if (is.null(fit_ptreat) & is.null(ptreat)) {
      fit_ptreat <- fit1
    }

    if (is.null(pred_ptreat) & is.null(ptreat)) {
      pred_ptreat <- pred1
    }

    # Generate formulas for first and second stage models
    formula1 <- .make_formula(outcome, c(grouping, basefeats))
    formula2 <- .make_formula(treatment, c("risk__", grouping, controls))


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

    train_ind <- data$fold__ == "train"
    ctl_train_ind <- train_ind & data[[treatment]] == 0
    trt_train_ind <- train_ind & data[[treatment]] == 1

    # Fit first-stage models
    if (is.null(resp_trt)) {
      m1_trt <- fit1(formula1, data[trt_train_ind, ], ...)

      if (calibrate) {
        data$resp_trt_pre_calib__ <- pred1(m1_trt, data, formula1)
        calib_formula <- .make_formula(outcome, "resp_trt_pre_calib__")

        data$resp_trt__ <-
          stats::predict.glm(stats::glm(calib_formula,
                                        data = data[trt_train_ind,],
                                        family = "binomial"),
                             data,
                             type = "response")
      } else {
        data$resp_trt__ <- pred1(m1_trt, data, formula1)
      }
    } else if (length(resp_trt) == nrow(data) | length(resp_trt) == 1) {
      m1_trt <- "Custom values of resp_trt provided"
      data$resp_trt__ <- resp_trt
    } else {
      stop("Bad specification of argument: resp_trt")
    }

    if (is.null(resp_ctl)) {
      m1_ctl <- fit1(formula1, data[ctl_train_ind, ], ...)

      if (calibrate) {
        data$resp_ctl_pre_calib__ <- pred1(m1_ctl, data, formula1)
        calib_formula <- .make_formula(outcome, "resp_ctl_pre_calib__")

        data$resp_ctl__ <-
          stats::predict.glm(stats::glm(calib_formula,
                                        data = data[ctl_train_ind,],
                                        family = "binomial"),
                             data,
                             type = "response")
      } else {
        data$resp_ctl__ <- pred1(m1_ctl, data, formula1)
      }
    } else if (length(resp_ctl) == nrow(data) | length(resp_ctl) == 1) {
      m1_ctl <- "Custom values of resp_ctl provided"
      data$resp_ctl__ <- resp_ctl
    } else {
      stop("Bad specification of argument: resp_ctl")
    }

    if (is.null(ptreat)) {
      m_ptrt <- fit_ptreat(formula, data[train_ind, ], ...)
      data$ptrt__ <- pred1(m_ptrt, data, formula)
    } else if (length(ptreat) == nrow(data) | length(ptreat) == 1) {
      # TODO(jongbin): Provide warning for cases when fit/pred_ptreat is
      # specified but ignored
      fit_ptreat <- "Overrided with custom values for ptreat"
      pred_ptreat <- "Overrided with custom values for ptreat"
      data$ptrt__ <- ptreat
    } else {
      stop("Bad specification of argument: ptreat")
    }

    data$risk__ <- logit(data[[paste0(risk, "__")]])

    ret <- list(data = data,
                m1_ctl = m1_ctl,
                m1_trt = m1_trt,
                risk_col = risk,
                treatment = treatment,
                outcome = outcome,
                grouping = grouping,
                features = basefeats,
                controls = controls,
                fit1 = fit1,
                pred1 = pred1,
                fit_ptreat = fit_ptreat,
                pred_ptreat = pred_ptreat,
                fit2 = fit2)

    class(ret) <- c("policy", class(ret))

    return(ret)
  }

#' Run test for unjustified disparate impact
#'
#' @param pol object of class policy
#' @param controls character vector of additional controls to consider in the
#'   second-stage model
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#'
#' @return tidy data frame of rad coefficients
#'
#' @export
compute_rad <-
  function(pol,
           controls = NULL,
           base_group = NULL,
           minority_groups = NULL) {
    # Input validation
    if (!("policy" %in% class(pol))) {
      stop("Expected object of class policy")
    }

    if (length(base_group) > 1) {
      stop("Specify a single base group.\n\tGot: ", base_group)
    }

    d <- pol$data
    group_col <- d[[pol$grouping]]
    members <- unique(group_col)

    check_groups <- sapply(c(base_group, minority_groups),
                           function(x) x %in% members)
    if (!all(check_groups)) {
      stop(sprintf("%s - not members of %s",
                   paste0(c(base_group, minority_groups)[!check_groups],
                          collapse = ","),
                   pol$grouping))
    }

    if (is.null(base_group)) {
      if (is.factor(group_col)) {
        base_group <- levels(group_col)[1]
      } else {
        base_group <- unique(group_col)[1]
      }
    }

    if (is.null(minority_groups)) {
      if (is.factor(group_col)) {
        minority_groups <- levels(group_col)[-1]
      } else {
        minority_groups <- unique(group_col)[-1]
      }
    }

    formula2 <- .make_formula(pol$treatment, c("risk__", pol$grouping, controls))
    test_df <- d[d$fold__ == "test", ]

    ret <- purrr::map_dfr(minority_groups, function(comp) {
      target_group_ind <- test_df[[pol$grouping]] %in% c(base_group, comp)

      tmp_df <- test_df[target_group_ind, ]
      tmp_df[[pol$grouping]] <- forcats::fct_drop(tmp_df[[pol$grouping]])
      tmp_df[[pol$grouping]] <- forcats::fct_relevel(tmp_df[[pol$grouping]],
                                                   base_group)

      coefs <- .pull_coefs(tmp_df, pol$treatment, pol$grouping,
                           c("risk__", controls),
                           fun = pol$fit2)

      coefs[grepl(pol$grouping, coefs$term), ]
    })

    return(ret)
  }
