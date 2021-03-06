#' Calculates the intersections of the tangents of a set of 
#' abscissae for a given density
#' @param abscissae vector of points at which to calculate tangents & intersections
#' @param density closure of the density of interest
#' @return y values of the intersection points of the tangents
calculate_tangents <- function(abscissae, density) {
  assertthat::assert_that(
    is.vector(abscissae),
    msg = "Abscissae is not a list"
  )
  assertthat::assert_that(
    typeof(density) == "closure",
    msg = "Density is not a closure"
  )
  
  tangents <- sapply(1:(length(abscissae) - 1), function(i) {
    
    return(
      (density(abscissae[i + 1]) - density(abscissae[i]) - abscissae[i + 1] * 
              numDeriv::grad(density, abscissae[i + 1], method = "simple") + 
         abscissae[i] * numDeriv::grad(density, abscissae[i], method = "simple")) / 
             (numDeriv::grad(density, abscissae[i], method = "simple") - 
                numDeriv::grad(density, abscissae[i + 1], method = "simple")))
  })
  
  return(sort(tangents))
}

#' Validates that a given density is log-concave for a set of points
#' for each subsequent pair of points, the function checks that 
#' the derivative conditions for log-concavity are met, ie the derivatives
#' are monotonically decreasing across the density's domain
#' @param abscissae sorted vector of points at which to evaluate the density
#' @param density closure of the density of interest
#' @return boolean indicating whether the function is log-concave
check_concavity <- function(abscissae, density) {
  assertthat::assert_that(
    is.vector(abscissae),
    msg = "Abscissae is not a list"
  )
  assertthat::assert_that(
    typeof(density) == "closure",
    msg = "Density is not a closure"
  )
  
  derivatives <- round(numDeriv::grad(density, abscissae), 6)
  
  pairwise_concavity <-
    sapply(1:(length(derivatives) - 1), function(i) {
      if (derivatives[i + 1] <= derivatives[i]) {
        return(TRUE)
      }
      
      return(FALSE)
      
    })
  
  # checks that concavity conditions are met for every pair
  if (sum(unlist(pairwise_concavity)) == (length(derivatives) - 1)) {
    return(TRUE)
  }
  
  return(FALSE)
}

#' Generate a sample from a given density using distr library
#' @param hull density of interest (a bounding hull in this context)
#' @return sample from the density
sample_from_hull <- function(hull, n_samples = 1) {
  assertthat::assert_that(
    typeof(hull) == "closure",
    msg = "Hull is not a closure"
  )
  
  assertthat::assert_that(
    is.numeric(n_samples) & n_samples > 0,
    msg = "Invalid n_samples parameter"
  )
  
  dist <- distr::AbscontDistribution(d = hull)
  rdist <- distr::r(dist)
  sample <- rdist(n_samples)
  return(sample)
}