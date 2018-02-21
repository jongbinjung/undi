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
