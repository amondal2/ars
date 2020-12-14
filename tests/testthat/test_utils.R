tol <- 1e-1
p_value_min <- .01
n <- 100

cat("\n")
test_that("test numerical derivative function", {
  print("Testing numerical derivative function from numDeriv library")
  # the deriv of the log normal dist d/dx log(e^(-x^2/2)/sqrt(2 π))= -x
  x <- 3
  result <-
    numDeriv::grad(get_log_density(dnorm), x, method = "simple")
  expected <- -x
  expect_equal(result, expected, tolerance = 1e-1)
  print("Test passed for derivative of log-normal density at x=3")
  
  # uniform should be 0
  x <- 0.3
  result <-
    numDeriv::grad(get_log_density(dunif), x, method = "simple")
  expected <- 0
  expect_equal(result, expected, tolerance = 1e-1)
  print("Test passed for derivative of uniform density at x=0.3")
})

cat("\n")
test_that("test distr sampling", {
  print("Testing numerical sampling function from distr library")
  # test normal case
  sample <- sample_from_hull(dnorm,n)
  p_val <- shapiro.test(sample)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for normal case")
  print("Shapiro-Wilks test of normality passed for distr sample")
})

cat("\n")
test_that("test generate_initial_abscissae", {
  # mode of a normal distribution is 0, with var 1
  print("Testing abscissae generating function")
  k <- 4
  mode <- 0
  scale <- 1
  expected_abscissae <-
    seq(mode - 2 * scale, mode + 2 * scale, length.out = k)
  expect_equal(expected_abscissae, generate_initial_abscissae(dnorm))
  print("Test passed for generate_initial_absisscae")
})