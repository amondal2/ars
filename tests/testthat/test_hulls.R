tol <- 1e-1

cat("\n")
testthat::test_that("test tangent calculations", {
  print("Testing tangent calculating function")
  # d/dx of log(std normal) = -x
  abscis <- c(-1, 0)
  
  #tangent to 0 is  y = -1/2 log(2 π)
  #and tangent line to -1 is  = x + 1 + 1/2 (-1 - log(2 π))
  # and they intersect at x=-.5
  expected <- -0.5
  result <- calculate_tangents(abscis, get_log_density(dnorm))
  expect_equal(result, expected, tolerance = tol)
  print("Test passed - tangent for log normal with starting abscissae at [-1, 0] found at x=-0.5")
  
})

cat("\n")
testthat::test_that("test concavity checking function", {
  print("Testing concavity checks for various distributions")
  #rnorm should pass
  expect(check_concavity(seq(-10, 10),
                         get_log_density(dnorm)),
         failure_message = "dnorm should pass check_concavity")
  
  print("Test passed - standard normal is log-concave")
  #rexp should pass
  expect(check_concavity(seq(.1, 10), get_log_density(dexp)),
         "dexp should pass check log concavity")
  print("Test passed - exponential is log-concave for rate=1")
  
  #cauchy should fail
  expect_false(check_concavity(seq(-10, 10), get_log_density(dcauchy)))
  print("Test passed - Cauchy is not log-concave for location=0, scale=1")
})