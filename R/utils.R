#' Calculates the derivative of the log density, using the chain rule
#' where f(x) = log(x) and g(x) = pdf, so d/dx f(g(x)) = g_prime(x)/g(x) 
#' @param density closure of the density of interest
#' @param x value of point at which to integrate
#' @param epsilon change value used to calculate derivative
#' @return value of derivative of log-density at a given point
deriv <- function(density, x, epsilon=.001) {
  g_prime <- (density(x + epsilon) - density(x)) / epsilon
  return(g_prime/density(x))
}

#' Numerically calculates the CDF for a given density
#' @param density closure of the density of interest
#' @return closure of the cdf which can be evaluated at a point
get_cdf <- function(density) {
  return (function(x) {
    # add lower bound of domain, calc number of points
    lower <- -25
    range <- seq(lower, x, length.out = 1000)
    return(sum(density(range))*((x-lower)/length(range)))
  })
}

#' Numerically inverts a given CDF 
#' @param prob probability value at which to invert the CDF
#' @param cdf closer of the cdf of interest
#' @return value of the inverse of the cdf for a given probability
eval_inverse_cdf <- function(prob, cdf) {
  assertthat::assert_that(prob >= 0 & prob <= 1)
  #TODO: use a better domain
  lower <- -100
  upper <- 100
  centered_cdf = function(x) {
    return(cdf(x) - prob)
  }
  return(uniroot(
    centered_cdf,
    lower = lower,
    upper = upper
  )$root)
}

#' Generate a sample from a given density using inverse CDF method
#' @param hull density of interest (a bounding hull in this context)
#' @return sample from the density
sample_from_hull <- function(hull) {
  return(eval_inverse_cdf(runif(1), get_cdf(hull)))
}