#' Entrypoint for 'ars' package. Runs an adaptative 
#' rejection sampling algorithm on a given probability density
#' and returns the specified number of samples.
#' @param density A function closure of the density to be sampled from
#' @param n_samples The number of samples to sample from the distribution
#' @param k Number of initial points from which to construct tangents
#' @param location starting point at which to initialize sampling points, default 0
#' @param scale value at which to space initial points default 1
#' @return samples numeric vector of samples from the specified distribution
#' @export
ars <- function(density, n_samples, k = 4, location = 0, scale = 1) {
  assertthat::assert_that(
    typeof(density) == "closure", 
    msg = "Density is not a function."
  )
  
  assertthat::assert_that(
    is.numeric(n_samples) & n_samples > 0,
    msg = "Invalid n_samples parameter"
  )
  
  
  log_density <- get_log_density(density)
  abscissae <- generate_initial_abscissae(density, location, scale, k)
 
  samples <- rep(0, n_samples)
  
  is_concave <- check_concavity(abscissae, log_density)
  assertthat::assert_that(
    is_concave == TRUE,
    msg = "Density is not log-concave for given set of points."
  )
  
  tangents <- calculate_tangents(abscissae, log_density)
  deriv_at_absc <- numDeriv::grad(log_density, abscissae, method="simple")
  
  upper_hull <- function(x) {
    intervals = findInterval(x, tangents)
    x_js = abscissae[intervals+1]
    x_js[is.na(x_js)] <- abscissae[length(abscissae)]
    derivs = deriv_at_absc[match(x_js, abscissae)]
    result = log_density(x_js) + (x - x_js) * derivs
    
    # use a vectorized version of the function for interop with R's integrate
    return(result)
  }
  
  exp_upper_hull <- function(x) {
    result <- exp(upper_hull(x))
    # zero out invalid densities
    result[density(x) <= 0] <- 0
    return(result)
  }
  
  integration_factor <- integrate(exp_upper_hull, -Inf, Inf)$value
  normalized_upper_hull <- function(x) {
    return(exp_upper_hull(x) / integration_factor)
  }
  
  lower_hull <- function(x) {
    intervals = findInterval(x, abscissae)
    x_js = c(-Inf,abscissae)[intervals+1]
    x_js_plus_1 = abscissae[intervals+1]
    res = (((x_js_plus_1 - x) * log_density(x_js) + 
          (x - x_js) * log_density(x_js_plus_1)) / (x_js_plus_1 - x_js))
    res[x < abscissae[1] | x > abscissae[length(abscissae)] | is.na(res)] <- -Inf
    return(res)
  }

  num_sampled <- 1
  
  while (num_sampled <= n_samples) {
    
    # we expect to be able to take larger batches as the hull
    # converges, which happens as we take more samples
    batch_size <- ceiling(num_sampled/2)
    
    sample <- sample_from_hull(normalized_upper_hull, batch_size)
    w <- runif(batch_size)
    
    # evaluate the hulls + density at the sampled point
    lower_bound <- lower_hull(sample)
    upper_bound <- upper_hull(sample)
    sample_log_density <- log_density(sample)
    
    to_accept <- sample[w<=exp(lower_bound - upper_bound)]
    num_to_accept <- length(to_accept)
    if (num_to_accept==batch_size) {
      # accept the entire sample
      samples[num_sampled:(num_sampled+num_to_accept-1)] <- to_accept
      num_sampled <- num_sampled + num_to_accept
    } else {
      # update tangents and points
      abscissae <- sort(c(abscissae, sample))
      abscissae <- abscissae[density(abscissae) > 0]

      # update all cached values that change when the abscissae change
      deriv_at_absc = numDeriv::grad(log_density, abscissae, method="simple")
      tangents <- calculate_tangents(abscissae, log_density)
      integration_factor <- integrate(exp_upper_hull, -Inf, Inf)$value
      
      # accept samples
      to_accept <- sample[w<=exp(sample_log_density - upper_bound)]
      num_to_accept <- length(to_accept)
      if(num_to_accept>0) {
        samples[num_sampled:(num_sampled+num_to_accept-1)] <- to_accept
        num_sampled <- num_sampled + num_to_accept
      }

    }
  }
  
  # check concavity before returning final results
  is_concave <- check_concavity(abscissae, log_density)
  assertthat::assert_that(
    is_concave == TRUE,
    msg = "Density is not log-concave for given set of points."
  )
  
  samples <- samples[1:n_samples] # trim extra samples from final batch
  return(samples)
}
