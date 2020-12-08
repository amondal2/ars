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
sample_from_hull <- function(hull) {
  dist <- distr::AbscontDistribution(d=hull) 
  rdist <- distr::r(dist)   
  sample <- rdist(1)
  return(sample)
}