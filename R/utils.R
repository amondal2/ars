#' Evaluates a function closure of the log of a given density function
#' @param density density of interest
#' @return closure representing the log density
get_log_density <- function(density) {
  assertthat::assert_that(
    typeof(density) == "closure",
    msg = "Density is not a function."
  )
  return(function(x) {
    return(log(density(x)))
  })
}

#' Generate the initial starting points for the algorithm
#' by calculating the mode of the function and then choosing
#' linearly spaced points around the mode
#' @param density density of interest (a bounding hull in this context)
#' @param location starting point of the mode-finding algorithm, default 0
#' @param scale width of the spacing of the points around the model, default 1
#' @param k number of points to initialize, default 4
#' @return vector of points
generate_initial_abscissae <-
  function(density,
           location = 0,
           scale = 1,
           k = 4) {
    assertthat::assert_that(
      typeof(density) == "closure",
      msg = "Density is not a function."
    )
    
    assertthat::assert_that(
      is.numeric(k) & k > 0,
      msg = "Invalid number of starting points"
    )
    
    assertthat::assert_that(
      is.numeric(location),
      msg = "Invalid location parameter"
    )
    
    assertthat::assert_that(
      is.numeric(scale) & scale > 0,
      msg = "Invalid scale parameter"
    )
    
    # nlm prints warnings if the density is 0 or undefined for x < 0 
    # (eg. exponential distribution), but
    # it doesn't affect the optimization so we suppress the warnings here
    mode <-
      suppressWarnings(nlm(function(x) {
        -1 * density(x)
      }, location)$estimate)
    abscissae <- seq(mode - 2 * scale, mode + 2 * scale, length.out = k)
    abscissae <- abscissae[density(abscissae) > 0]
    return(abscissae)
  }