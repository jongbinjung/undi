
<!-- README.md is generated from README.Rmd. Please edit that file -->
undi
====

[![Travis-CI Build Status](https://travis-ci.org/jongbinjung/undi.svg?branch=master)](https://travis-ci.org/jongbinjung/undi) [![Coverage Status](https://img.shields.io/codecov/c/github/jongbinjung/undi/master.svg)](https://codecov.io/github/jongbinjung/undi?branch=master)

Installation
------------

You can install undi from github with:

``` r
# install.packages("devtools")
devtools::install_github("jongbinjung/undi")
```

Example
-------

``` r
library(undi)
library(tidyverse)
inv_logit <- stats::binomial()$linkinv

N <- 2000
set.seed(1)
data <- tibble(id = rep(1:N)) %>%
  mutate(
  m = rnorm(N, 0, 1),
  x = factor(rep(c("red", "blue"), N / 2)),
  x = fct_relevel(x, "red"),
  z = rnorm(N, m, 2),
  c = rnorm(N, m, 2),
  e1 = rnorm(N, 0, 0.05),
  e2 = rnorm(N, 0, 0.05),
  risk = 1 + 0.05 * (x == "blue") + 0.25 * z + c + e1,
  a = inv_logit(risk + e2) > .5,
  y = inv_logit(risk) > .5
  )
```

``` r
example_policy <-
  policy(a ~ x + z + c, 
         data, 
         outcome = "y", 
         cv.folds = 2, 
         distribution = "bernoulli")

# Initial RAD computation
compute_rad(example_policy)
#>    term estimate std.error statistic p.value controls
#> 1 xblue   0.3236    0.2392     1.353  0.1763   risk__

# Sensitivity with uniform parameters
sensitivity(example_policy, 
            q = .45, 
            dp = log(1.8), 
            d0 = log(2), 
            d1 = log(1.5))
#>    term estimate std.error statistic p.value controls
#> 3 xblue  0.02014  0.004904     4.106 4.1e-05   risk__

# Sensitivity with parameters assigned to levels of the grouping variable
sensitivity(example_policy, 
            q = c(.45, .55), 
            dp = c(-log(1.2), log(1.5)), 
            d0 = c(-log(3), log(3)), 
            d1 = c(0, log(1.8))
  )
#>    term estimate std.error statistic p.value controls
#> 3 xblue  0.01237  0.005585     2.215 0.02681   risk__
```
