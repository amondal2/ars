#' Evaluates a function closure of the log of a given density function
#' @param density density of interest 
#' @return closure representing the log density
get_log_density <- function(density) {
  return(function(x) {
    return(log(density(x)))
  })
}

#' Generate a sample from a given density using distr library
#' @param hull density of interest (a bounding hull in this context)
#' @return sample from the density
sample_from_hull <- function(hull, n_samples=1) {
  dist <- distr::AbscontDistribution(d=hull) 
  rdist <- distr::r(dist)   
  sample <- rdist(n_samples)
  return(sample)
}

#' Generate the initial starting points for the algorithm
#' by calculating the mode of the function and then choosing
#' linearly spaced points around the mode
#' @param density density of interest (a bounding hull in this context)
#' @param location starting point of the mode-finding algorithm, default 0
#' @param scale width of the spacing of the points around the model, default 1
#' @param k number of points to initialize, default 4
#' @return vector of points
generate_initial_abscissae <- function(density, location=0, scale=1, k=4) {
  log_density <- get_log_density(density)
  
  # nlm prints warnings if the density is 0 or undefined for x < 0, but
  # it doesn't affect the optimization so we suppress the warnings here
  mode <- suppressWarnings(nlm(function(x) {-1*log_density(x)}, location)$estimate)
  abscissae <- seq(mode-2*scale, mode+2*scale, length.out = k)
  abscissae <- abscissae[density(abscissae) > 0]
  return(abscissae)
}