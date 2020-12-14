tol = 1e-1
p_value_min = .01
n = 100


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