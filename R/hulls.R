#' Calculates the intersections of the tangents of a set of 
#' abscissae for a given density
#' @param abscissae vector of points at which to calculate tangents & intersections
#' @return y values of the intersection points of the tangents
calculate_tangents <- function(abscissae, density) {
  tangents <- sapply(1:(length(abscissae) - 1), function(i) {
    
    return((log(density(abscissae[i + 1])) - log(density(abscissae[i])) - abscissae[i + 1] * 
              deriv(density, abscissae[i + 1]) + abscissae[i] * deriv(density, abscissae[i])) / 
             (deriv(density, abscissae[i]) - deriv(density, abscissae[i + 1])))
  })
  
  return(tangents)
}