#' Calculates the intersections of the tangents of a set of 
#' abscissae for a given density
#' @param abscissae vector of points at which to calculate tangents & intersections
#' @param density closure of the density of interest
#' @return y values of the intersection points of the tangents
calculate_tangents <- function(abscissae, density) {
  tangents <- sapply(1:(length(abscissae) - 1), function(i) {
    
    return((log(density(abscissae[i + 1])) - log(density(abscissae[i])) - abscissae[i + 1] * 
              deriv(density, abscissae[i + 1]) + abscissae[i] * deriv(density, abscissae[i])) / 
             (deriv(density, abscissae[i]) - deriv(density, abscissae[i + 1])))
  })
  
  return(tangents)
}

#' Validates that a given density is log-concave for a set of points
#' for each subsequent pair of points, the function checks that 
#' the derivative conditions for log-concavity are met
#' @param abscissae vector of points at which to evaluate the density
#' @param density closure of the density of interest
#' @return boolean indicating whether the function is log-concave
check_concavity <- function(abscissae, density) {
  pairwise_concavity <- sapply(1:(length(abscissae) - 1), function(i) {
    deriv_x_i <- deriv(density, abscissae[i])
    deriv_x_i1 <- deriv(density,abscissae[i+1])
    # todo this is sometimes NaN
    
    if(deriv_x_i * deriv_x_i1 > 0) {
      if(deriv_x_i - deriv_x_i1 > 0) {
        return(TRUE)
      }
    } else if (deriv_x_i * deriv_x_i1 < 0) {
      if(deriv_x_i > 0 & deriv_x_i1 < 0) {
        return(TRUE)
      }
      return(FALSE)
      
    }
  })
  
  # checks that concavity conditions are met for every pair
  if(sum(unlist(pairwise_concavity)) == (length(abscissae)-1)) {
    return(TRUE)
  }
  
  return(FALSE)
}