.make_formula <- function(y, vars) {
  stats::as.formula(paste(y, "~", paste(vars, collapse = "+")))
}

logit <- stats::binomial()$linkfun
inv_logit <- stats::binomial()$linkinv


#' Overwrites a list with named arguments
#'
#' Returns a list that is the intersection of
#' \code{defaults} and \code{...} where named
#' values in \code{...} overwrite values in \code{defaults}.
#'
#' One or more of the elements in \code{...} can be a
#' \code{pairlist} (the structure in which dotted arguments
#' are stored) and this list will get expanded and substituted
#' into \code{defaults} too. This means that the same named argument
#' could be listed more than twice, in which case the right-most value
#' is used. See below for details.
#'
#' @param defaults list of default values.
#'   Argument always be specified by name not position
#' @param ... a comma separated set of named
#'
#' \code{overwrite_list(defaults = list(a = 1, b = '2'), c = 1, a = '3')} will
#' return list(a = '1', b = '2', c = 1)
#'
#' To demonstrate the operation when one of the arguments is a \code{pairlist},
#' note that both of the functions below produce the same result. You can either pass the
#' elipsis (...) directly, or you can extract the arguments from the \code{call}
#' and pass them later.
#'
#' \code{
#'   f = function(...) {
#'     overwrite_list(defaults = list(a = 1, b = 'b'), b = 2, ...)
#'   }
#' }
#'
#' \code{
#'   f = function(...) {
#'     # Extract passed arguments from call
#'     dot_args = match.call(expand.dots = F)$...
#'     overwrite_list(defaults = list(a = 1, b = 'b'), b = 2, dot_args)
#'   }
#' }
#'
#' \code{f(k = 5)} will return \code{list(a = 1, b = 2, k = 5)}.
#'
#' \code{f(b = 10)} will return \code{list(a = 1, b = 10)}. Note
#' how the named argument \code{b} appears 3 times in the \code{overwrite_list}
#' call; once in the defaults, once explicitly specified for overwriting (\code{b = 2}),
#' and once indirectly specified via the arguments to \code{f} (\code{b = 10}).
#'
#' @export
overwrite_list <- function(defaults, ...) {
  dot_args <- list(...)

  # Expand pairlists (ie lists of arguments)
  dot_args =
    lapply(seq_along(dot_args), function(i) {
      if (is.pairlist(dot_args[[i]])) {
        do.call(c,dot_args[i])
      } else {
        dot_args[i]
      }}) %>%
    unlist(recursive = F)

  for (v in seq_along(dot_args)) {
    defaults[[names(dot_args)[v]]] <- dot_args[[v]]
  }

  defaults
}


#' Given a policy object, compute and return the appropriate risk column
.get_risk_col <- function(pol) {
  # Input validation
  if (!("policy" %in% class(pol))) {
    stop("Expected object of class policy")
  }

  logit(pol$data[[paste0(pol$risk_col, "__")]])
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

#' Extract sensitivity parameters from a list of parameters for optimization
#'
#' @param params vector of parameters
#' @param free_params (Optional) logical vector of length 8 indicating which
#'   parameters are specified in params. \code{length(params)} should equal
#'   \code{sum(free_params == TRUE)}.
#' @param fixed_param_values (Optional) vector specifing the values of the fixed
#'   parameters (ie when \code{free_params == FALSE})
#' @param q_range if true, the second parameter defines the log odds ratio
#'   between q for base and minority groups
#' @param allow_sgv logical; whether to allow for subgroup validity; i.e., if
#'   \code{TRUE}, the delta parameters (\code{dp}, \code{d0}, \code{d1}) will be
#'   allowed to vary between base/minority groups, but if \code{FALSE}, a single
#'   value for each delta parameter will be used for each base/minority pair
#'
#' @return named list of parameters
.extract_params <- function(params,
                            free_params = rep(T, 8),
                            fixed_param_values = NULL,
                            q_range = FALSE,
                            allow_sgv = FALSE) {

  all_params = numeric(8)
  all_params[free_params] = params
  all_params[!free_params] = fixed_param_values

  if (q_range) {
    all_params[2] = inv_logit(logit(all_params[1]) + all_params[2])
  }

  if (!allow_sgv) {
    all_params[4] <- all_params[3]
    all_params[6] <- all_params[5]
    all_params[8] <- all_params[7]
  }

  list(qb  = all_params[1],
       qm  = all_params[2],
       ab  = all_params[3],
       am  = all_params[4],
       d0b = all_params[5],
       d0m = all_params[6],
       d1b = all_params[7],
       d1m = all_params[8])
}

#' Compute AUC given predictions and corresponding labels
#'
#' @param pred scalar predictions
#' @param label binary labels
#' @param ret_num logical, whether to return a numerical value (TRUE) or a
#'   formatted string (FALSE: default)
#'
#' @return scalar AUC
.compute_auc <- function(pred, label, ret_num = FALSE) {
  if (length(unique(label)) == 2) {
    p <- ROCR::prediction(pred, label)
    auc <- ROCR::performance(p, "auc")

    if (ret_num) {
      return(unlist(auc@y.values))
    } else {
      return(format(unlist(auc@y.values) * 100))
    }
  } else {
    if (ret_num) {
      return(-1)
    } else {
      return("-")
    }
  }
}

#' Down-sample a data set
#'
#' @param d original data frame
#' @param n parameter for down-sampling; if between 0 and 1, treated as a
#'   fraction, if larger than 1, rounded and treated as the number of rows to
#'   sample
#' @param verbose whether to print down-sampling message
#'
#' @return down-sampled data frame
.down_sample <- function(d, n, verbose = TRUE) {
  N <- nrow(d)

  if (n == 1) {
    return(d)
  }

  if (n <= 0) {
    stop("Down sampling parameter must be greater than 0.\n\tGot: ", n)
  }

  if (n < 1) {
    ret <- dplyr::sample_frac(d, size = n)
  } else {
    ret <- dplyr::sample_n(d, size = n)
  }

  if (verbose) {
    message("Down-sampled data to ", nrow(ret), "/", N,
            " (", format(nrow(ret)/N * 100), "%) rows")
  }

  return(ret)
}

#' Safely scale/center covariates, avoiding NaNs
#'
#' @param x a numeric matrix (like object)
#'
#' @return scaled matrix with attributes \code{scaled:center}
#'   and \code{scaled:scale} indicating the transformation performed.
.safe_scale <- function(x) {

  means = colMeans(x, na.rm = T)
  sds = apply(x, 2, function(c) stats::sd(c, na.rm = T))

  badcols = sds < .Machine$double.eps*2

  # Don't scale constant columns
  sds[badcols] = 1
  means[badcols] = 0

  scale(x, center = means, scale = sds)
}


#' Extract group members from a column of groups
#'
#' @param group_col column of groups
#'
#' @return \code{levels(group_col)} if \code{group_col} is a factor, and
#'   \code{unique(group_col)} otherwise
.get_groups <- function(group_col) {
  if (is.factor(group_col)) {
    return(levels(group_col))
  } else {
    return(unique(group_col))
  }
}
