#' Entrypoint for 'ars' package. Runs an adaptative
#' rejection sampling algorithm on a given probability density
#' and returns the specified number of samples.
#' @param density A function closure of the density to be sampled from.
#' Input functions must have similar behavior to dnorm, dexp, etc.
#' @param n_samples The number of samples to sample from the distribution
#' @param location starting point at which to initialize sampling points, 
#' default 0. Should be chosen so that `location` is near the mode 
#' of the distribution.
#' @param scale value at which to space initial points (default 1). 
#' Should be chosen such that mode - (2*scale), mode+(2*scale) is within
#'  the support of the probability density.
#' @return samples numeric vector of samples from the specified distribution
#' @export
ars <-
  function(density,
           n_samples = 10,
           location = 0,
           scale = 1) {
    
    num_starting_abscissae <- 4
    
    assertthat::assert_that(
      typeof(density) == "closure",
      msg = "Density is not a function."
    )
    
    assertthat::assert_that(
      is.numeric(n_samples) & n_samples > 0,
      msg = "Invalid n_samples parameter"
    )
    
    log_density <- get_log_density(density)
    
    abscissae <- generate_initial_abscissae(
      density, 
      location, 
      scale, 
      num_starting_abscissae
    )
    samples <- rep(0, n_samples)
    
    is_concave <- check_concavity(abscissae, log_density)
    assertthat::assert_that(
      is_concave == TRUE,
      msg = "Density is not log-concave for given set of points."
    )
    
    tangents <- calculate_tangents(abscissae, log_density)
    deriv_at_absc <- numDeriv::grad(
      log_density, 
      abscissae, 
      method = "simple"
    )
    
    upper_hull <- function(x) {
      intervals <- findInterval(x, tangents)
      # vectorized calculation across all abscissae
      x_j <- abscissae[intervals + 1]
      x_j[is.na(x_j)] <- abscissae[length(abscissae)]
      derivs <- deriv_at_absc[match(x_j, abscissae)]
      result <- log_density(x_j) + (x - x_j) * derivs
      
      return(result)
    }
    
    exp_upper_hull <- function(x) {
      result <- exp(upper_hull(x))
      # zero out invalid densities
      result[density(x) <= 0] <- 0
      return(result)
    }
    
    # cache the normalizing constant to not recalculate 
    # for every sample
    integration_factor <- integrate(exp_upper_hull, -Inf, Inf)$value
    normalized_upper_hull <- function(x) {
      return(exp_upper_hull(x) / integration_factor)
    }
    
    lower_hull <- function(x) {
      # vectorized calculation of lower hull across abscissae
      intervals <- findInterval(x, abscissae)
      x_j <- c(-Inf, abscissae)[intervals + 1]
      x_j_plus_1 <- abscissae[intervals + 1]
      result <- (((x_j_plus_1 - x) * log_density(x_j) +
                (x - x_j) * log_density(x_j_plus_1)
      ) / (x_j_plus_1 - x_j))
      result[x < abscissae[1] |
          x > abscissae[length(abscissae)] | 
            is.na(result)
      ] <- -Inf
      return(result)
    }
    
    num_sampled <- 1
    
    while (num_sampled <= n_samples) {
      # we expect to be able to take larger batches as the hull
      # converges, which happens as we take more samples
      batch_size <- ceiling(num_sampled / 2)
      
      sample <- sample_from_hull(normalized_upper_hull, batch_size)
      w <- runif(batch_size)
      
      # evaluate the hulls + density at the sampled point
      lower_bound <- lower_hull(sample)
      upper_bound <- upper_hull(sample)
      sample_log_density <- log_density(sample)
      
      to_accept_mask <- (w <= exp(lower_bound - upper_bound))
      to_accept <- sample[to_accept_mask]
      num_to_accept <- length(to_accept)
      if (num_to_accept == batch_size) {
        # accept the entire sample
        samples[num_sampled:(num_sampled + num_to_accept - 1)] <-
          to_accept
        num_sampled <- num_sampled + num_to_accept
      } else {
        # update tangents and points
        abscissae <- sort(c(abscissae, sample[!to_accept_mask]))
        abscissae <- abscissae[density(abscissae) > 0]
        
        # update all cached values that change when the abscissae change
        deriv_at_absc <- numDeriv::grad(log_density, abscissae, method = "simple")
        tangents <- calculate_tangents(abscissae, log_density)
        integration_factor <-
          integrate(exp_upper_hull, -Inf, Inf)$value
        
        # accept samples
        to_accept <- sample[w <= exp(sample_log_density - upper_bound)]
        num_to_accept <- length(to_accept)
        if (num_to_accept > 0) {
          samples[num_sampled:(num_sampled + num_to_accept - 1)] <- to_accept
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
    
    # trim extra samples from final batch
    samples <- samples[1:n_samples] 
    return(samples)
  }
