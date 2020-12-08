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
  
  log_density <- get_log_density(density)
  abscissae <- seq(-5, 5, length.out = k)
 
  samples <- rep(0, n_samples)
  
  if (density(-Inf) == 0) {
    # todo optimize how x1 is selected
    abscissae[1] <- -5
  }
  
  if (density(Inf) == 0) {
    # todo optimize how x_k is selected
    abscissae[k] <- 5
  }
  
  #is_concave <- check_concavity(abscissae, log_density)
  # assertthat::assert_that(
  #   is_concave == TRUE,
  #   msg = "Density is not log-concave for given set of points."
  # )
  
  tangents <- calculate_tangents(abscissae, log_density)
  
  upper_hull <- function(x) {
    x_j <- abscissae[k]
    # TODO fix tangent error
    for (i in 1:length(tangents)) {
      if (tangents[i] > x[1]) {
        x_j <- abscissae[i]
        break
      }
    }
    return(log_density(x_j) + (x - x_j) * numDeriv::grad(log_density, x_j))
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
        (abscissae[j + 1] - x) * log_density(abscissae[j]) + 
          (x - abscissae[j]) * log_density(abscissae[j + 1])) / (abscissae[j + 1] - abscissae[j]))
  }
  
  num_sampled <- 1
  while (num_sampled <= n_samples) {
    sample <- sample_from_hull(normalized_upper_hull)
    w <- runif(1)
    
    # evaluate the hulls + density at the sampled point
    lower_bound <- lower_hull(sample)
    upper_bound <- upper_hull(sample)
    sample_log_density <- log_density(sample)
    
    if (w <= exp(lower_bound - upper_bound)) {
      # accept the sample
      samples[num_sampled] <- sample
      num_sampled <- num_sampled + 1
    } else {
      # update tangents and points, check concavity conditions
      # are still met
      abscissae <- sort(c(abscissae, sample))
      #is_concave <- check_concavity(abscissae, log_density)
      
      # assertthat::assert_that(
      #   is_concave == TRUE,
      #   msg = "Density is not log-concave for given set of points."
      # )
      
      tangents <- calculate_tangents(abscissae, log_density)
      
      if (w <= exp(sample_log_density - upper_bound)) {
        # accept the sample only if this condition is met
        samples[num_sampled] <- sample
        num_sampled <- num_sampled + 1
      }
    }
  }
  return(samples)
}
