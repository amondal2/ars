library(assertthat)

deriv <- function(density, x, epsilon=.001) {
  #calculates the derivative of the log density, using the chain rule
  # where f(x) = log(x) and g(x) = pdf, so
  # d/dx f(g(x)) = g_prime(x)/g(x) 
  g_prime = (density(x + epsilon) - density(x)) / epsilon
  return(g_prime/density(x))
}

get_cdf <- function(density) {
  return (function(x) {
    # add lower bound of domain, calc number of points
    lower <- -25
    range <- seq(lower, x, length.out = 1000)
    return(sum(density(range))*((x-lower)/length(range)))
  })
}

eval_inverse_cdf <- function(prob, cdf) {
  # given a probability and a cdf function,
  # find the x value that corresponds to the given prob
  
  #assert_that(0<=prob)
  #assert_that(prob<=1)
  
  #TODO: use a better domain
  lower <- -100
  upper <- 100
  centered_cdf = function(x) {
    return(cdf(x)-prob)
  }
  return(uniroot(centered_cdf, lower=lower, upper=upper)$root)
}


ars <- function(density, n_samples, k=5) {
  # todo validation: 1. check density is log-concave, N-samples > 0
  x_k <- seq(-25, 25, length.out = k)

  if(density(-Inf) == 0) {
    # todo optimize how x1 is selected
    x_k[1] <- -25
  }

  if(density(Inf) == 0) {
    # todo optimize how x_k is selected
    x_k[k] <- 25
  }

  z_j <- sapply(1:(length(x_k) - 1), function(i) {

    # check if this should be symmetric

    return((log(density(x_k[i + 1])) - log(density(x_k[i])) - x_k[i + 1] * 
              deriv(density, x_k[i + 1]) + x_k[i] * deriv(density, x_k[i])) / 
                (deriv(density, x_k[i]) - deriv(density, x_k[i + 1])))
  })
  
  


  u_k <- function(x) {
    x_j <- x_k[k]
    for(i in 1:length(z_j)) {
      if(z_j[i] > x) {
        x_j <- x_k[i]
        break
      }
    }
    return(log(density(x_j)) + (x-x_j)*deriv(density, x_j))
  }
  
  exp_u_k <- function(x) {
    return(exp(u_k(x)))
  }
  
  s_k <- function(x) {
    #todo figure out bounds
    return(exp(sapply(x, u_k))/integrate(exp_u_k, -Inf, Inf)$value)
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
  
  x_star <- eval_inverse_cdf(runif(1), get_cdf(s_k))
  w <- runif(1)
  if (w <= exp(l_k(x_star) - u_k(x_star))) {
    print("accept")
  } else {
    #todo figure out if we need derivative
    if (w <= exp(density(x_star) - u_k(x_star))) {
      print("accept")
    } else {
      print("reject")
    }
    
  }
  return(x_star)

}


