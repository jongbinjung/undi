#' Find parameter values that min/maximize sensitivity results via optim
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
#' @param range_q_ratio (Optional) 2D vector. If set, the minority and base
#'   values of q will not be allowed to vary independently, but instead will be
#'   constrained to vary by the given range of log odds. ie q_minority =
#'   inv.logit(logit(q_base) + u), where u is in range_q_ratio
#' @param allow_sgv logical; whether to allow for subgroup validity; i.e., if
#'   \code{TRUE}, the delta parameters (\code{dp}, \code{d0}, \code{d1}) will be
#'   allowed to vary between base/minority groups, but if \code{FALSE}, a single
#'   value for each delta parameter will be used for each base/minority pair
#' @param controls vector of legitimate controls to use; the ones specified
#'   within the policy object will be used if not specified
#' @param optim_control list of control parameters passed to \code{optim}
#' @param include_benchmark logical; whether to include the two extreme
#'   benchmark test results (default: FALSE)
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
#' @inheritParams rad_control
#' @export
optimsens <-
  function(pol,
           range_q = c(0, 1),
           range_dp = c(0, log(2)),
           range_d0 = c(0, log(2)),
           range_d1 = c(0, log(2)),
           base_group = NULL,
           minority_groups = NULL,
           range_q_ratio = NULL,
           allow_sgv = FALSE,
           controls = NULL,
           fit_fn = "logit_coef",
           optim_control = list(),
           include_benchmark = FALSE,
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
    params_upper[2] <- range_q_ratio[2]
  }

  free_params <- abs(params_upper - params_lower) > 2*.Machine$double.eps

  if (!allow_sgv) {
    # Treat minority delta parameters as "fixed" values
    # TODO: Might want to eliminate "magic number" indices
    free_params[c(4, 6, 8)] <- FALSE
  }

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

  # Create grid of min/max option and minority groups for iteration
  iter_grid <- list(minor = minority_groups,
                    sgn = c(-1, 1)) %>%
    purrr::cross_df() %>%
    distinct()

  # Generate initial values
  # TODO(jongbin): Allow the user to specify initial values?
  params_init <- furrr::future_map_dfr(1:nrow(iter_grid),
                                       function(ip) {
    it <- iter_grid[ip, ]

    tag_ <- paste(it$minor, ifelse(it$sgn > 0, "min", "max"), sep = "_")

    init_ <- stats::runif(length(params_lower),
                          min = params_lower,
                          max = params_upper)

    tibble(tag = tag_,  params = list(init_))
  })

  # Optimize for min/max over each minority group
  optim_res <- furrr::future_map_dfr(1:nrow(iter_grid),
                                     function(ip) {
    it <- iter_grid[ip, ]
    minor <- it$minor
    sgn <- it$sgn
    tag_ <- paste(minor, ifelse(sgn > 0, "min", "max"), sep = "_")
    pars <- params_init %>%
      filter(tag == tag_) %>%
      pull("params") %>%
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
        allow_sgv = allow_sgv,
        controls = controls,
        fit_fn = fit_fn,
        naive_se = FALSE,
        verbose = verbose
      ),
      lower = params_lower[free_params],
      upper = params_upper[free_params],
      method = 'L-BFGS-B',
      control = optim_control
    )

    # Convert parameters back to 8-D vector
    opt$par = unlist(.extract_params(
      opt$par, free_params,
      pars[!free_params],
      !is.null(range_q_ratio)))

    tibble(minor = minor,
           tag = tag_,
           optim = list(opt))
  })

  # Extract final results from optimized parameters
  coefs <-
    furrr::future_map_dfr(1:nrow(optim_res), function(ip) {
      minor <- optim_res[ip, ][["minor"]]
      tag_ <- optim_res[ip, ][["tag"]]
      params <- optim_res[ip, ] %>%
        mutate(pars = purrr::map(optim, "par")) %>%
        pull("pars") %>%
        `[[`(1)

      fn <- .get_optim_fn(pol, sgn = sgn,
                          compare = c(base_group, minor),
                          controls = controls, naive_se = TRUE,
                          return_scalar = FALSE)
      ret <- fn(params)

      ret <- ret %>%
        mutate(tag = tag_,
               pars = list(params))

      ret
    })

  base_case <- bind_rows(
    lapply(
      minority_groups,
      function(group) {
        sensitivity(pol, 0, 0, 0, 0, compare = c(base_group, group),
                    controls = controls, fit_fn = fit_fn)
      }))

  base_case$method <- "rad"

  if (include_benchmark) {
    base_bm <-  bind_rows(
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

    base_case <- bind_rows(base_bm, base_case)
  }

  ret <- list(
    results = coefs,
    optim = optim_res,
    base_case = base_case,
    base_group = base_group)

  class(ret) <- c("optimsens", "sens", class(ret))

  return(ret)
}


#' Find parameter values that min/maximize sensitivity results via grid search
#'
#' Within specified values of sensitivity parameters, find the ones that achieve
#' minimum/maximum sensitivity results
#'
#' @param pol object of class \code{policy}
#' @param qs vector q values to search
#' @param dps values to search for change in log-odds of treat = 1 if u = 1
#' @param d0s values to search for change in log-odds of response = 1 if treat =
#'   0 and u = 1
#' @param d1s values to search for change in log-odds of response = 1 if treat =
#'   1 and u = 1
#' @param params_grid (Optional) data frame with columns
#'   \code{qb, qm, ab, am, d0b, d0m, d1b, d1m} indicating the parameter combinations
#'   to be searched over. If specified, the parameter range argments (\code{qs} through \code{d1s})
#'   are ignored.
#' @param base_group (Optional) single group that acts as the pivot/base; by
#'   default, if the grouping variable is a factor, set to the first level,
#'   otherwise set to the first of sorted unique values
#' @param minority_groups (Optional) groups to compare to the base group; by
#'   default, set to every unique value other than the base group
#' @param allow_sgv logical; whether to allow for subgroup validity; i.e., if
#'   \code{TRUE}, the delta parameters (\code{dp}, \code{d0}, \code{d1}) will be
#'   allowed to vary between base/minority groups, but if \code{FALSE}, a single
#'   value for each delta parameter will be used for each base/minority pair
#' @param controls vector of legitimate controls to use; the ones specified
#'   within the policy object will be used if not specified
#' @param include_benchmark logical; whether to include the two extreme
#'   benchmark test results (default: FALSE)
#' @param verbose whether or not to print debug messages (0 = none, 1 = results
#'   only, 2 = everything)
#'
#' @return a list-type object of class \code{gridsens} with the following
#'   elements \item{results}{\code{tidy} dataframe of second-stage model
#'   coefficients after searching for min/max values across specified sensitivy
#'   parameter values, independently for each minority group}
#'   \item{grid}{results from full grid}
#'   \item{base_case}{result from \code{compute_rad} on base policy with
#'   specified groups and controls} \item{base_group}{base group used in
#'   analysis}
#'
#' @inheritParams rad_control
#' @export
gridsens <-
  function(pol,
           qs = c(0, .4, .5, .6, 1),
           dps = c(0, log(2)),
           d0s = c(0, log(2)),
           d1s = c(0, log(2)),
           params_grid = NULL,
           base_group = NULL,
           minority_groups = NULL,
           allow_sgv = FALSE,
           controls = NULL,
           fit_fn = "logit_coef",
           include_benchmark = FALSE,
           verbose = TRUE) {
  # Input validation
  if (!("policy" %in% class(pol))) {
    stop("Expected object of class policy")
  }

  if (length(base_group) > 1) {
    stop("Specify a single base group.\n\tGot: ", base_group)
  }

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

  if (is.null(params_grid)) {
    params_grid <- .get_params_grid(qs,
                                    dps,
                                    d0s,
                                    d1s,
                                    minority_groups,
                                    allow_sgv)
  } else {
    # Check params_grid for validity
    if (is.null(params_grid$qb) ||
        any((params_grid$qb < 0) | (params_grid$qb > 1)) ||
        is.null(params_grid$qm) ||
        any((params_grid$qm < 0) | (params_grid$qm > 1)) ||
        is.null(params_grid$ab) ||
        any(params_grid$ab < 0) ||
        is.null(params_grid$am) ||
        any(params_grid$am < 0) ||
        is.null(params_grid$d0b) ||
        any(params_grid$d0b < 0) ||
        is.null(params_grid$d0m) ||
        any(params_grid$d0m < 0) ||
        is.null(params_grid$d1b) ||
        any(params_grid$d1b < 0) ||
        is.null(params_grid$d1m) ||
        any(params_grid$d1m < 0)) {
      stop('params_grid is misspecified')
    }
  }

  # Optimize for min/max over each minority group
  grid_res <- furrr::future_map_dfr(1:nrow(params_grid),
                                    .progress = verbose,
                                    function(ip) {
      params <- params_grid[ip, ]

      if (verbose >= 2) {
        cat(sprintf("grid: %s\n", .format_pars(params)))
      }

      ret <- sensitivity(
        pol,
        compare = c(base_group, params$minority),
        q = c(params$qb, params$qm),
        dp = c(params$ab, params$am),
        d0 = c(params$d0b, params$d0m),
        d1 = c(params$d1b, params$d1m),
        fit_fn = fit_fn,
        controls = controls,
        naive_se = FALSE,
        verbose = FALSE
        ) %>%
        mutate(pars = list(params), minor = params$minority)
    })

  grid_opt <- grid_res %>%
    group_by(term, controls) %>%
    mutate(max = max(estimate), min = min(estimate)) %>%
    filter(estimate == max | estimate == min) %>%
    # Remove duplicate maxima/minima (by selecting top-most result)
    group_by(estimate, add = TRUE) %>%
    slice(1) %>%
    mutate(tag = paste(minor, ifelse(estimate == max, "max", "min"), sep = "_"))

  # Extract final results from optimized parameters
  coefs <- furrr::future_map_dfr(1:nrow(grid_opt), function(ip) {
    params <- grid_opt[ip, ] %>%
      pull("pars") %>%
      `[[`(1)

    minor <- grid_opt[ip, ]$minor
    tag <- grid_opt[ip, ]$tag

    ret <- sensitivity(
      pol,
      compare = c(base_group, minor),
      q = c(params$qb, params$qm),
      dp = c(params$ab, params$am),
      d0 = c(params$d0b, params$d0m),
      d1 = c(params$d1b, params$d1m),
      controls = controls,
      naive_se = TRUE,
      fit_fn = fit_fn,
      verbose = FALSE
      ) %>%
      mutate(pars = list(params))

    ret$tag <- tag
    ret
  })

  base_case <- bind_rows(
    lapply(
      minority_groups,
      function(group) {
        sensitivity(pol, 0, 0, 0, 0, compare = c(base_group, group),
                    fit_fn = fit_fn,
                    controls = controls)
      }))

  base_case$method <- "rad"

  if (include_benchmark) {
    base_bm <-  bind_rows(
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

    base_case <- bind_rows(base_bm, base_case)
  }

  ret <- list(
    results = coefs,
    grid = grid_res,
    base_case = base_case,
    base_group = base_group)

  class(ret) <- c("gridsens", "sens", class(ret))

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
#' @param tag string to tag output with (usefull for parallel output)
#' @param controls vector of controls
#' @param q_range if true, the second parameter defines the log odds ratio
#'   between q for base and minority
#' @param allow_sgv logical; whether to allow for subgroup validity; i.e., if
#'   \code{TRUE}, the delta parameters (\code{dp}, \code{d0}, \code{d1}) will be
#'   allowed to vary between base/minority groups, but if \code{FALSE}, a single
#'   value for each delta parameter will be used for each base/minority pair
#' @param naive_se whether or not to compute "naive" std.errors in sensitivity;
#'   FALSE by default, to avoid unnecessary computation, but should be computed
#'   for final results once extreme values have been identified
#' @param verbose whether or not to print debug messages (0 = none, 1 = results
#'   only, 2 = everything)
#' @param return_scalar logical, whether to return a single scalar values (TRUE)
#'   or to return the full result from \code{sensitivity}
#'
#' @inheritParams rad_control
#' @return Function that will return coefficient on minority group
.get_optim_fn <-
  function (pol,
            sgn,
            compare,
            free_params = NULL,
            fixed_param_values = NULL,
            tag = "fit",
            controls = NULL,
            q_range = FALSE,
            allow_sgv = FALSE,
            naive_se = FALSE,
            fit_fn = "logit_coef",
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
                         q_range = q_range,
                         allow_sgv = allow_sgv)
    qb  <- p$qb
    qm  <- p$qm
    ab  <- p$ab
    am  <- p$am
    d0b <- p$d0b
    d0m <- p$d0m
    d1b <- p$d1b
    d1m <- p$d1m

    if (verbose >= 2) {
      cat(sprintf("%s: %s\n", tag, .format_pars(p)))
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
      fit_fn = fit_fn,
      verbose = FALSE
    )

    if (return_scalar) {
      ret <- ret[["estimate"]]

      if (verbose) {
        cat(sprintf("\t%s coef: %.4f\n", tag, ret))
      }

      ret <- ret * sgn
    }

    ret
  }
}


#' Generate data frame grid of parameter values, given possible values
#'
#' @param qs vector q values to search
#' @param dps values to search for change in log-odds of treat = 1 if u = 1
#' @param d0s values to search for change in log-odds of response = 1 if treat =
#'   0 and u = 1
#' @param d1s values to search for change in log-odds of response = 1 if treat =
#'   1 and u = 1
#' @param allow_sgv logical; whether to allow for subgroup validity; i.e., if
#'   \code{TRUE}, the delta parameters (\code{dp}, \code{d0}, \code{d1}) will be
#'   allowed to vary between base/minority groups, but if \code{FALSE}, a single
#'   value for each delta parameter will be used for each base/minority pair
#'
#' @return data frame where each row represents a unique combination of possible
#'   parameter values
.get_params_grid <- function(qs, dps, d0s, d1s, mg, allow_sgv) {
  if (allow_sgv) {
    params_list <-
      list(
        minority = mg,
        qb = qs,
        qm = qs,
        ab = dps,
        am = dps,
        d0b = d0s,
        d0m = d0s,
        d1b = d1s,
        d1m = d1s
      )
    params_grid <- purrr::cross_df(params_list)
  } else {
    params_list <-
      list(
        minority = mg,
        qb = qs,
        qm = qs,
        dp = dps,
        d0 = d0s,
        d1 = d1s
      )
    params_grid <- purrr::cross_df(params_list) %>%
      mutate(
        ab = dp,
        am = dp,
        d0b = d0,
        d0m = d0,
        d1b = d1,
        d1m = d1
      ) %>%
      select(-dp,-d0,-d1)
  }

  # TODO(jongbin): Better way of filtering "insignificant" parameter combos?
  params_grid <- params_grid %>%
    mutate(
      maxparam = pmax(ab, am, d0b, d0m, d1b, d1m),
      minq = pmin(qb, qm),
      maxq = pmax(qb, qm)
    ) %>%
    filter(maxparam != 0,
           minq != 1,
           maxq != 0) %>%
    select(-maxparam,-minq,-maxq) %>%
    distinct()

  return(params_grid)
}
