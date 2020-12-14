context("ars")

tol = 1e-1
p_value_min = .01
n = 100

test_that("test input validation", {
  density <- 5
  expect_error(ars(density, 5))
  
  n_samples <- "test"
  expect_error(ars(dnorm, n_samples))
})

test_that("test resulting distribution", {
  sample = ars(dnorm,n)
  p_val = shapiro.test(sample)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for normal case")
  
})


test_that("test exponential distribution", {
  x <- rexp(n)
  y <- ars(dexp, n, location=1)
  p_val <- ks.test(x,y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for exponential case")
})

test_that("test chisq distribution", {
  df <- 5
  x <- rchisq(n,df)
  y <- ars(function(x){dchisq(x, df)}, n, location=1)
  p_val <- ks.test(x,y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for chisq case")
})

test_that("test gamma distribution", {
  alpha <- 6
  beta <- 9
  x <- rgamma(n,alpha,beta)
  y <- ars(function(x){dgamma(x, alpha,beta)}, n, location=1)
  p_val <- ks.test(x,y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for gamma case")
})

test_that("test beta distribution", {
  alpha <- 6
  beta <- 9
  x <- rbeta(n,alpha,beta)
  y <- ars(function(x){dbeta(x, alpha,beta)}, n, location=.5, scale=.1)
  p_val <- ks.test(x,y)$p.value
  expect(p_val >= p_value_min, "rnorm p-value below limit for gamma case")
})






