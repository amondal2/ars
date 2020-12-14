tol = 1e-1
p_value_min = .01
n = 100


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
  n
  sample = sample_from_hull(dnorm,n)
  p_val = shapiro.test(sample)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for normal case")
})

test_that("test generate_initial_abscissae", {
  # mode of a normal distribution is 0, with var 1
  k <- 4
  mode <- 0
  scale <- 1
  expected_abscissae <- seq(mode-2*scale, mode+2*scale, length.out = k)
  expect_equal(expected_abscissae, generate_initial_abscissae(dnorm))
})