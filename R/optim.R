#' Find parameter values that min/maximize sensitivity results
#'
#' Within specified range of sensitivity parameters, find the ones that achieve
#' minimum/maximum sensitivity results
#'
#' @param pol object of class \code{policy}
#' @param range_q 2D vector specifying min/max value of p(u = 1 | x)
#' @param range_dp 2D vector specifying min/max value of change in log-odds of
#'   treat = 1 if u = 1
#' @param range_d0 2D vector specifying min/max value of change in log-odds of
#'   response = 1 if treat = 0 and u = 1
#' @param range_d1 2D vector specifying min/max value of change in log-odds of
#'   response = 1 of treat = 1 and u = 1
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#' @param range_q_ratio (Optional) 2D vector. If set, the minority and majority
#'   values of q will not be allowed to vary independently, but instead will be
#'   constrained to vary by the given range of log odds. ie q_minority =
#'   inv.logit(logit(q_majority) + u), where u is in range_q_ratio
#' @param controls vector of legitimate controls to use; the ones specified
#'   within the policy object will be used if not specified
#' @param include_benchmark logical; whether to include the two extreme
#'   benchmark test results (default: FALSE)
#' @param transform_params logical; whether to transform parameters such with
#'   logit(q) and log(delta) --- might improve optimization conditioning
#' @param verbose whether or not to print debug messages (0 = none, 1 = results
#'   only, 2 = everything)
#' @param debug logical flag, if TRUE, returns a list of results and the
#'   expanded data frame used to fit model
#'
#' @details If any of the \code{range_} arguments are set to a single value
#'   (instead of a 2D vector), the corresponding paramter will be fixed and not
#'   explored for the optimization
#'
#' @return a list-type object of class \code{optimsens} with the following
#'   elements \item{results}{\code{tidy} dataframe of second-stage model
#'   coefficients after searching for min/max values across specified sensitivy
#'   parameter ranges, independently for each minority group}
#'   \item{optim}{nested data frame where the \code{$optim} column contains the
#'   optimization results} \item{base_case}{result from \code{compute_rad} on
#'   base policy with specified groups and controls} \item{base_group}{base
#'   group used in analysis}
#'
#' @export
optimsens <-
  function(pol,
           range_q = c(-Inf, Inf),
           range_dp = c(-Inf, log(log(2))),
           range_d0 = c(-Inf, log(log(2))),
           range_d1 = c(-Inf, log(log(2))),
           base_group = NULL,
           minority_groups = NULL,
           range_q_ratio = NULL,
           controls = NULL,
           include_benchmark = FALSE,
           transform_params = TRUE,
           verbose = TRUE,
           debug = FALSE) {
  # Input validation
  if (!("policy" %in% class(pol))) {
    stop("Expected object of class policy")
  }

  if (length(base_group) > 1) {
    stop("Specify a single base group.\n\tGot: ", base_group)
  }

  if (!is.null(range_q_ratio) && (length(range_q_ratio) != 2)) {
    stop("range_q_ratio should be NULL or a vector of length 2")
  }

  # range_q = c(-Inf, Inf)
  # range_dp = c(-Inf, log(log(2)))
  # range_d0 = c(-Inf, log(log(2)))
  # range_d1 = c(-Inf, log(log(2)))
  # range_q = c(0, 1)
  # range_dp = c(0, log(2))
  # range_d0 = c(0, log(2))
  # range_d1 = c(0, log(2))

  group_col <- pol$data[[pol$grouping]]
  members <- unique(group_col)

  check_groups <- sapply(c(base_group, minority_groups),
                         function(x) x %in% members)
  if (!all(check_groups)) {
    stop(sprintf("%s - not members of %s",
                 paste0(c(base_group, minority_groups)[!check_groups],
                        collapse = ","),
                 pol$grouping))
  }

  # Optimization hyper parameters
  params_lower <- c(
    min(range_q),   # Min P(u = 1 | base )
    # Min P(u = 1 | minority ) or
    # Min log-odds(u = 1 | base) - log-odds(u = 1 | minority)
    min(range_q),
    min(range_dp),  # Min alpha_base
    min(range_dp),  # Min alpha_minority
    min(range_d0),  # Min delta0_base
    min(range_d0),  # Min delta0_minority
    min(range_d1),  # Min delta1_base
    min(range_d1)   # Min delta1_minority
  )

  if (!is.null(range_q_ratio)) {
    params_lower[2] = range_q_ratio[1]
  }

  params_upper <- c(
    max(range_q),   # Max P(u = 1 | base )
    # Max P(u = 1 | minority ) or
    # Max log-odds(u = 1 | base) - log-odds(u = 1 | minority)
    max(range_q),
    max(range_dp),  # Max alpha_base
    max(range_dp),  # Max alpha_minority
    max(range_d0),  # Max delta0_base
    max(range_d0),  # Max delta0_minority
    max(range_d1),  # Max delta1_base
    max(range_d1)   # Max delta1_minority
  )

  if (!is.null(range_q_ratio)) {
    params_upper[2] = range_q_ratio[2]
  }

  free_params = abs(params_upper - params_lower) > 2*.Machine$double.eps

  if (any(params_upper[free_params] < params_lower[free_params])) {
    stop("Upper value in range must be equal to or larger than lower value")
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

  # Generate initial values
  # TODO(jongbin): Allow the user to specify initial values?
  params_init <- foreach(minor = minority_groups,
                         .combine = dplyr::bind_rows) %:%
    foreach(sgn = c(-1, 1), .combine = dplyr::bind_rows) %dopar% {
      tag_ <- paste(minor, ifelse(sgn > 0, "min", "max"), sep = "_")

      init_ <- stats::runif(length(params_lower),
                            min = pmax(params_lower, -100),
                            max = pmin(params_upper, 0))

      dplyr::tibble(tag = tag_,  params = list(init_))
  }

  # Optimize for min/max over each minority group
  optim_res <- foreach(minor = minority_groups,
                       .combine = dplyr::bind_rows) %:%
    foreach(sgn = c(-1, 1), .combine = dplyr::bind_rows) %dopar% {
      tag_ <- paste(minor, ifelse(sgn > 0, "min", "max"), sep = "_")
      pars <- params_init %>%
        dplyr::filter(tag == tag_) %>%
        dplyr::pull("params") %>%
        `[[`(1)

      opt = stats::optim(
        pars[free_params],
        .get_optim_fn(
          pol,
          sgn = sgn,
          free_params = free_params,
          fixed_param_values = pars[!free_params],
          compare = c(base_group, minor),
          tag = tag_,
          q_range = !is.null(range_q_ratio),
          controls = controls,
          naive_se = FALSE,
          transform_params = transform_params,
          verbose = verbose
        ),
        lower = params_lower[free_params],
        upper = params_upper[free_params],
        method = 'L-BFGS-B'
      )

      # Convert parameters back to 8-D vector
      opt$par = unlist(.extract_params(
        opt$par, free_params,
        pars[!free_params],
        !is.null(range_q_ratio)))

      dplyr::tibble(minor = minor,
                    tag = tag_,
                    optim = list(opt))
    }

  # Extract final results from optimized parameters
  coefs <-
    foreach(ip = 1:nrow(optim_res), .combine = dplyr::bind_rows) %dopar% {
    minor <- optim_res[ip, ][["minor"]]
    tag_ <- optim_res[ip, ][["tag"]]
    params <- optim_res[ip, ] %>%
      dplyr::mutate(pars = purrr::map(optim, "par")) %>%
      dplyr::pull("pars") %>%
      `[[`(1)

    fn <- .get_optim_fn(
      pol,
      sgn = sgn,
      compare = c(base_group, minor),
      controls = controls,
      naive_se = TRUE,
      transform_params = transform_params,
      return_scalar = FALSE
    )

    ret <- fn(params)
    ret$tag <- tag_
    ret
    }

  base_case <- compute_rad(
    pol,
    controls = controls,
    base_group = base_group,
    minority_groups = minority_groups
  )
  base_case$method <- "rad"

  if (include_benchmark) {
    base_bm <-  dplyr::bind_rows(
      compute_bm(pol,
                 base_group = base_group,
                 minority_groups = minority_groups),
      compute_bm(
        pol,
        base_group = base_group,
        minority_groups = minority_groups,
        kitchen_sink = TRUE
      )
    )
    base_bm$method <- "bm"

    base_case <- dplyr::bind_rows(base_bm, base_case)
  }

  ret <- list(
    results = coefs,
    optim = optim_res,
    base_case = base_case,
    base_group = base_group)

  class(ret) <- c("optimsens", class(ret))

  return(ret)
}


#' Generate subroutine for computing sensitized race coefficients for single
#' minority group v. whites
#'
#' @param pol policy object
#' @param sgn sign integer to multiply on scalar return value; used for
#'   controlling max/min optimization
#' @param compare vector of length 2, specifying the two groups to compare
#' @param free_params a logical vector indicating which of the (8) input
#'   parameters are free to vary. If NULL, all parameters are allowed to vary
#' @param fixed_param_values A vector defining the values of the fixed params.
#'   \code{length(fixed_param_values)} should equal \code{sum(free_params == F)}
#' @param transform_params logical; whether to transform parameters such with
#'   logit(q) and log(delta) --- might improve optimization conditioning
#' @param tag string to tag output with (usefull for parallel output)
#' @param controls vector of controls
#' @param q_range if true, the second parameter defines the log odds ratio
#'   between q for majority and minority
#' @param naive_se whether or not to compute "naive" std.errors in sensitivity;
#'   FALSE by default, to avoid unnecessary computation, but should be computed
#'   for final results once extreme values have been identified
#' @param verbose whether or not to print debug messages (0 = none, 1 = results
#'   only, 2 = everything)
#' @param return_scalar logical, whether to return a single scalar values (TRUE)
#'   or to return the full result from \code{sensitivity}
#'
#' @return Function that will return coefficient on minority group
.get_optim_fn <-
  function (pol,
            sgn,
            compare,
            free_params = NULL,
            fixed_param_values = NULL,
            transform_params = TRUE,
            tag = "fit",
            controls = NULL,
            q_range = FALSE,
            naive_se = FALSE,
            verbose = TRUE,
            return_scalar = TRUE) {

  # Validate input
  if (length(compare) != 2) {
    stop("Can only get optim fn for comparison of two groups, got ",
         length(compare))
  }

  if (is.null(free_params)) {
    free_params = rep(T, 8)
  }

  if (length(free_params) != 8 || !is.logical(free_params)) {
    stop("free_params must be a logical vector of length 8")
  }

  if (sum(!free_params) != length(fixed_param_values)) {
    stop(paste0("The length of fixed_param_values does not equal",
                "the number of fixed parameters indicated by free_params"))
  }

  function(params) {

    # Unpack parameters
    p <- .extract_params(params,
                         free_params = free_params,
                         fixed_param_values = fixed_param_values,
                         transform_params = transform_params,
                         q_range = q_range)

    qb  <- p$qb
    qm  <- p$qm
    ab  <- p$ab
    am  <- p$am
    d0b <- p$d0b
    d0m <- p$d0m
    d1b <- p$d1b
    d1m <- p$d1m

    if (verbose >= 2) {
      cat(sprintf(paste("%s: q=%.2f/%.2f, e(a)=%.2f/%.2f,",
                        "e(d0)=%.2f/%.2f, e(d1)=%.2f/%.2f\n"),
                  tag, qb, qm, exp(ab), exp(am), exp(d0b), exp(d0m),
                  exp(d1b), exp(d1m)))
    }

    ret <- sensitivity(
      pol,
      compare = compare,
      q = c(qb, qm),
      dp = c(ab, am),
      d0 = c(d0b, d0m),
      d1 = c(d1b, d1m),
      controls = controls,
      naive_se = naive_se,
      verbose = FALSE
    )

    if (return_scalar) {
      ret <- ret[["estimate"]]

      if (verbose) {
        cat(sprintf("  %s coef: %.4f\n", tag, ret))
      }

      ret <- ret * sgn
    }

    ret
  }
}

