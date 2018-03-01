context("run")


# Generate top-level test data --------------------------------------------
N <- 2000
set.seed(1)
data <- data.frame(id = rep(1:N)) %>%
  dplyr::mutate(
    m = rnorm(N, 0, 1),
    x = factor(rep(c("red", "blue"), N / 2)),
    z = rnorm(N, m, 2),
    c = rnorm(N, m, 2),
    e1 = rnorm(N, 0, 0.05),
    e2 = rnorm(N, 0, 0.05),
    risk = 1 + 0.05 * (x == "blue") + 0.25 * z + c + e1,
    a = inv_logit(risk + e2) > .5,
    y = inv_logit(risk) > .5
  )

# Test basic constructor functionality ------------------------------------
test_that("undi can generate a trivial object with custom values", {
  undi_obj <- undi(a ~ x + z + c, data, outcome = "y",
                   ptreat = .5,
                   # Dummy function to skip model fitting
                   fit2 = function(f, d) lm(1 ~ 1, d),
                   resp_ctl = 0,
                   resp_trt = 1)

  expect_is(undi_obj, "undi")
  expect_equal(undi_obj$data$ptrt__, rep(.5, nrow(data)))
  expect_equal(undi_obj$data$resp_ctl__, rep(0, nrow(data)))
  expect_equal(undi_obj$data$resp_trt__, rep(1, nrow(data)))
})
