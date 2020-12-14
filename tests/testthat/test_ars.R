context("ars")

tol <- 1e-1
p_value_min <- .01
n <- 100

testthat::test_that("test input validation", {
  print("Testing ars input validation.")
  density <- 5
  expect_error(ars(density, 5))
  print("Density is not a function - test expect_error passed")
  
  n_samples <- "test"
  expect_error(ars(dnorm, n_samples))
  print("n_samples is not a number - test expect_error passed")
  
})

cat("\n")
testthat::test_that("test normal distribution", {
  print("Testing sample generation from standard normal")
  sample <- ars(dnorm,n)
  p_val <- shapiro.test(sample)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for normal case")
  print("Shapiro-Wilks test of normality passed for ars sample")
  
})

cat("\n")
test_that("test exponential distribution", {
  print("Testing sample generation from exponential distribution")
  x <- rexp(n)
  y <- ars(dexp, n, location=1)
  p_val <- ks.test(x,y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for exponential case")
  print("Kolmogorv-Smirnov test passed for ars sample vs sample from rexp")
  
})

cat("\n")
test_that("test chisq distribution", {
  print("Testing sample generation from chi-squared distribution (df=5)")
  df <- 5
  x <- rchisq(n, df)
  y <- ars(function(x) {
    dchisq(x, df)
  }, n, location = 1)
  p_val <- ks.test(x, y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for chisq case")
  print("Kolmogorv-Smirnov test passed for ars sample vs sample from rchisq")
  
})

cat("\n")
test_that("test gamma distribution", {
  print("Testing sample generation from gamma distribution (a=6, b=9)")
  alpha <- 6
  beta <- 9
  x <- rgamma(n, alpha, beta)
  y <- ars(function(x) {
    dgamma(x, alpha, beta)
  }, n, location = 1)
  p_val <- ks.test(x, y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for gamma case")
  print("Kolmogorv-Smirnov test passed for ars sample vs sample from rgamma")
})

cat("\n")
test_that("test beta distribution", {
  print("Testing sample generation from beta distribution (a=6, b=9)")
  alpha <- 6
  beta <- 9
  x <- rbeta(n, alpha, beta)
  y <-
    ars(function(x) {
      dbeta(x, alpha, beta)
    }, n, location = .5, scale = .1)
  p_val <- ks.test(x, y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for beta case")
  print("Kolmogorv-Smirnov test passed for ars sample vs sample from rbeta")
})
