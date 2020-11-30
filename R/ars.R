deriv <- function(density, x) {
  epsilon = .01
  return (log(density(x + epsilon) - log(density(x)))) / epsilon
}


ars <- function(density, n_samples, k=5) {
  # todo validation: 1. check density is log-concave, N-samples > 0
  x_k <- seq(-10, 10, length.out = 5)

  if(density(-Inf) == 0) {
    # todo optimize how x1 is selected
    x_k[1] <- -10
  }

  if(density(Inf) == 0) {
    # todo optimize how x_k is selected
    x_k[k] <- 10
  }

  z_j <- sapply(1:(length(x_k) - 1), function(i) {

    # check if this should be symmetric

    return(
      log(density(x_k[i + 1])) - log(density(x_k[i])) - x_k[i + 1] * deriv(density, x_k[i +
                                                                                          1]) + x_k[i] * deriv(density, x_k[i])
    ) / (deriv(density, x_k[i]) - deriv(density, x_k[i + 1]))
  })


  u_k <- function(x) {
    for(i in 1:length(z_j)) {
      if(z_j[i] > x) {
        x_j <- x_k[i]
        break
      }
    }
    return(log(density(x_j)) + (x-x_j)*deriv(density, x_j))
  }

  s_k <- function(x) {
    #todo figure out bounds
    return(exp(u_k(x))/integrate(u_k, z_j[1], z_j[k]))
  }

  l_k <- function(x) {
    #todo add check for x < x_1 or x > x_k
    for(i in 1:(length(x_k)-1)) {
      if(x_k[i+1]>x) {
        j <- i
        break
      }
    }
    return(
      ((x_k[j+1]-x)*log(density(x_k[j]))+(x-x_k[j])*log(density(x_k[j+1])))/(x_k[j+1]-x_k[j])
    )
  }


}


