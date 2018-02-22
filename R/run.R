#' Run test for unjustified disparate impact
#'
#' @param formula a formula in the form of \code{treatment ~ grouping_variable +
#'   other controls} where the LHS is the treatment column, first element on the
#'   RHS is the grouping variable (e.g., Race), and the remainder of the RHS
#'   specifies the form of controls
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
#' @param pred1 a function of the form f(model, data) used for generating
#'   predictions from the first-stage model; predictions should be on the scale
#'   that is appropriate for the second stage model, i.e., the second stage
#'   model will be fit with the formula \code{treatment ~ risk + controls}
#' @param fit2 a function of the form f(formula, data) used for fitting the
#'   second-stage model; using \code{glm} with \code{family = binomial} by
#'   default
#' @param seed random seed to use
#' @param ... additional arguments passed to first-stage model fitting function,
#'   \code{fit1}
#'
#' @details The first-stage model, which estimates \code{risk} is always trained
#'   on the subset of data where \code{treatment == 1}, i.e., \code{risk} is
#'   generally considered to be P(outcome = 1 | treatment = 1)
#'
#' @return undi object
#' @export
# TODO: Update documentation for return type undi
undi <-
  function(formula,
           data,
           outcome,
           controls = NULL,
           train = 0.5,
           fit1 = NULL,
           pred1 = NULL,
           fit2 = NULL,
           seed = 1234,  # TODO: set seed randomly
           ...) {
    set.seed(seed)

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
      pred1 <- function(m, d) gbm::predict.gbm(m, d,
                                               gbm::gbm.perf(m),
                                               type = "link")
    }

    if (is.null(fit2)) {
      fit2 <- function(f, d) stats::glm(f, d, family = stats::binomial)
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
        data$fold__ <- data[[train]]
      } else {
        stop("train should either be between 0 and 1, or column name")
      }
    } else if (length(train) == nrow(data)) {
      data$fold__ <- ifelse(train, "train", "test")
    } else {
      stop("Wrong specification of argument train; see ?undi")
    }

    train_df <- data[data$fold__ == "train", ]
    treated_train_ind <- train_df[[treatment]] == 1

    # Fit first-stage model
    m1 <- fit1(formula1, train_df[treated_train_ind, ], ...)
    data$risk__ <- pred1(m1, data)

    test_df <- data[data$fold__ == "test", ]

    coefs <- .pull_coefs(test_df, treatment, grouping,
                c("risk__", controls),
                fun = fit2) %>%
      dplyr::filter(grepl(grouping, term))

    ret <- list(data = data,
                m1 = m1,
                treatment = treatment,
                outcome = outcome,
                grouping = grouping,
                features = basefeats,
                controls = controls,
                fit1 = fit1,
                pred1 = pred1,
                fit2 = fit2,
                coefs = coefs)

    class(ret) <- c("undi", class(ret))

    return(ret)
  }
