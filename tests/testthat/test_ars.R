context("ars")

tol = 1e-1
p_value_min = .01

test_that("test input validation", {
  density <- 5
  expect_error(ars(density, 5))
  
  n_samples <- "test"
  expect_error(ars(dnorm, n_samples))
})

test_that("test resulting distribution", {
  n = 100
  sample = ars(dnorm,n)
  p_val = shapiro.test(sample)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for normal case")
  
})


test_that("test two samples came from the same distribution", {
  n <- 500
  x <- rexp(n)
  y <- ars(dexp, n, location=1)
  p_val <- ks.test(x,y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for normal case")
})


test_that("test deriv", {
  # the deriv of the log normal dist d/dx log(e^(-x^2/2)/sqrt(2 π))= -x
  x=3
  result = numDeriv::grad(get_log_density(dnorm), x, method="simple")
  expected = -x
  expect_equal(result,expected, tolerance=1e-1)
  
  #unif should be 0
  x=.3
  result = numDeriv::grad(get_log_density(dunif), x, method="simple")
  expected = 0
  expect_equal(result,expected, tolerance=1e-1)
})

test_that("test distr sampling", {
  # test normal case
  n = 1000
  sample = sample_from_hull(dnorm,n)
  p_val = shapiro.test(sample)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for normal case")
})

test_that("test calc tangents", {
  # d/dx of log(std normal) = -x
  abscis = c(-1,0)
  
  #tangent to 0 is  y = -1/2 log(2 π)
  #and tangent line to -1 is  = x + 1 + 1/2 (-1 - log(2 π)) 
  # and they intersect at x=-.5
  expected = -.5
  result = calculate_tangents(abscis, get_log_density(dnorm))
  expect_equal(result,expected, tolerance=tol)
  
})

test_that("check log concavity", {
  #rnorm should pass
  expect(check_concavity(seq(-10,10), get_log_density(dnorm)),
         failure_message = "dnorm should pass check_concavity")
  
  #rexp should pass
  expect(check_concavity(seq(.1,10), get_log_density(dexp)),
         "dexp should pass check log concavity")
  
  #cauchy should fail 
  expect_false(check_concavity(seq(-10,10), get_log_density(dcauchy)))
})

test_that("test generate_initial_abscissae", {
  # mode of a normal distribution is 0, with var 1
  k <- 4
  mode <- 0
  scale <- 1
  expected_abscissae <- seq(mode-2*scale, mode+2*scale, length.out = k)
  expect_equal(expected_abscissae, generate_initial_abscissae(dnorm))
})