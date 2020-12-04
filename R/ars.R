#' Entrypoint for 'ars' package. Runs an adaptative 
#' rejection sampling algorithm on a given probability density
#' and returns the specified number of samples.
#' @param density A function closure of the density to be sampled from
#' @param n_samples The number of samples to sample from the distribution
#' @param k Number of initial points from which to construct tangents
#' @return samples numeric vector of samples from the specified distribution
#' @export
ars <- function(density, n_samples, k = 5) {
  # should we be able to pass in params & domain?
  assertthat::assert_that(
    typeof(density) == "closure", 
    msg = "Density is not a function."
  )
  
  assertthat::assert_that(
    is.numeric(n_samples) & n_samples > 0,
    msg = "Invalid n_samples parameter"
  )
  
  abscissae <- seq(-25, 25, length.out = k)
  samples <- vector()
  
  if (density(-Inf) == 0) {
    # todo optimize how x1 is selected
    abscissae[1] <- -25
  }
  
  if (density(Inf) == 0) {
    # todo optimize how x_k is selected
    abscissae[k] <- 25
  }
  
  tangents <- calculate_tangents(abscissae, density)
  
  upper_hull <- function(x) {
    x_j <- abscissae[k]
    for (i in 1:length(tangents)) {
      if (tangents[i] > x[1]) {
        x_j <- abscissae[i]
        break
      }
    }
    return(log(density(x_j)) + (x - x_j) * deriv(density, x_j))
  }
  
  exp_upper_hull <- function(x) {
    return(exp(upper_hull(x)))
  }
  
  normalized_upper_hull <- function(x) {
    #todo figure out bounds
    return(exp(sapply(x, upper_hull)) / integrate(exp_upper_hull, -Inf, Inf)$value)
  }
  
  lower_hull <- function(x) {
    if(x < abscissae[1] | x > abscissae[k]) return(-Inf)
    
    for (i in 1:(length(abscissae) - 1)) {
      if (abscissae[i + 1] > x) {
        j <- i
        break
      }
    }
    return(
      (
        (abscissae[j + 1] - x) * log(density(abscissae[j])) + 
          (x - abscissae[j]) * log(density(abscissae[j + 1]))) / (abscissae[j + 1] - abscissae[j]))
  }
  
  while (length(samples) < n_samples) {
    sample <- sample_from_hull(normalized_upper_hull)
    w <- runif(1)
    
    # evaluate the hulls + density at the sampled point
    lower_bound <- lower_hull(sample)
    upper_bound <- upper_hull(sample)
    log_density <- log(density(sample))
    
    
    assertthat::assert_that(
      lower_bound <= log_density & upper_bound >= log_density,
      msg="Density is not log-concave."
    )
    
    if (w <= exp(lower_bound - upper_bound)) {
      # accept the sample
      samples <- sort(c(samples, sample))
    } else {
      # update tangents and points
      abscissae <- sort(c(abscissae, sample))
      tangents <- calculate_tangents(abscissae, density)
      
      if (w <= exp(log_density - upper_bound)) {
        # accept the sample only if this condition is met
        samples <- sort(c(samples, sample))
      }
    }
  }
  return(samples)
}
