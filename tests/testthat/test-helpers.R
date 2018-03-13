context("helpers")


# Test .make_formula() ----------------------------------------------------
test_that(".make_formula generates proper single-variable formulas", {
  target <- y ~ x
  generated <- .make_formula("y", "x")

  expect_equal(generated, target)
})

test_that(".make_formula generates proper multi-variable formulas", {
  target <- y ~ x1 + x2 + x3
  generated <- .make_formula("y", c("x1", "x2", "x3"))

  expect_equal(generated, target)
})

# Test .extract_features() ------------------------------------------------
test_that(".extract_features creates proper list from linear formula", {
  target <- list(treat = "y", group = "x", feats = c("a", "b", "c"))
  generated <- .extract_features(y ~ x + a + b + c)

  expect_equal(generated, target)
})

test_that(".extract_features extracts interactions", {
  target <- list(treat = "y", group = "x", feats = c("a", "b:z", "a:b:c"))
  generated <- .extract_features(y ~ x + a + b:z + a:b:c)

  expect_equal(generated, target)
})


# Test .expand_params() ---------------------------------------------------
test_that(".expand_params works for atomic parameters", {
  g <- factor(rep(c("a", "b", "c"), 3))
  p <- 1.23

  target <- rep(p, length(g))
  generated <- .expand_params(g, p)

  expect_equal(generated, target)
})

test_that(".expand_params works for parameters mapped to levels", {
  g <- factor(sample(c("a", "b", "c"), 10, replace = TRUE))
  p <- seq(1, length(levels(g)))

  target <- p[as.numeric(g)]
  generated <- .expand_params(g, p)

  expect_equal(generated, target)
})

test_that(".expand_params works for parameters with equal size as group", {
  g <- factor(sample(c("a", "b", "c"), 10, replace = TRUE))
  p <- rnorm(10)

  target <- p
  generated <- .expand_params(g, p)

  expect_equal(generated, target)
})

test_that(".expand_params throws errors for wrong specification", {
  g <- sample(c("a", "b", "c"), 10, replace = TRUE)

  # Level-length assignment for non-factor groupings
  p <- seq(1, length(unique(g)))
  expect_error(.expand_params(g, p))

  # Parameters of bad length
  p <- seq(1, length(g) + 1)
  expect_error(.expand_params(g, p))

  p <- seq(1, length(g) - 1)
  expect_error(.expand_params(g, p))

  p <- c(1, 2)
  expect_error(.expand_params(g, p))

  p <- NULL
  expect_error(.expand_params(g, p))

  p <- 1:1e4
  expect_error(.expand_params(g, p))
})

# Test .extract_params() --------------------------------------------------
test_that(".extract_params maps as expected", {
  params <- stats::runif(8)
  target <- list(qb  = params[1],
                 qm  = params[2],
                 ab  = params[3],
                 am  = params[4],
                 d0b = params[5],
                 d0m = params[6],
                 d1b = params[7],
                 d1m = params[8])
  generated <- .extract_params(params)

  expect_equal(generated, target)
})


# Test .compute_auc() -----------------------------------------------------
test_that(".compute_auc returns '-' for single-value labels", {
  preds <- runif(10)
  labels <- rep(TRUE, 10)
  generated <- .compute_auc(preds, labels)
  expect_equal(generated, "-")
})

test_that(".compute_auc returns .5 for single-value predictors", {
  preds <- rep(.5, 10)
  labels <- sample(c(TRUE, FALSE), 10, replace = TRUE)
  generated <- .compute_auc(preds, labels)
  expect_equal(generated, .5)
})
