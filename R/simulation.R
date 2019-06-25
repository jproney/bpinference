#' Simulate multi-type Markov branching process
#'
#' @param e_mat The matrix of birth events that can occur in the brnaching process. Dimensions \code{nevents} x \code{ntypes}
#' @param p_vec A vector containing the parent type for each of the birth events in \code{e_mat}. Dimensions \code{nevents} x 1
#' @param r_vec A vector containing the rate at which each of the birth events in \code{e_mat} occurs. Dimensions \code{nevents} x 1
#' @param z0_vec The initial population vector at time 0. Dimensions \code{ntypes} x 1
#' @param times The timepoints at which to record the population
#' @return A matrix with the population vectors at each timepoint. Dimensions \code{ntimes} x \code{ntypes}
#'
#' @export
bp <- function(e_mat, r_vec, p_vec, z0_vec, times)
{
  # z0_vec = inital population vector, times = vector of timepoints to record
  tf <- times[length(times)]
  ntypes <- length(z0_vec)
  nevents <- nrow(e_mat)
  zt_mat <- matrix(0, ncol = ntypes, nrow = length(times))  #store populations at different times
  t <- 0
  zt_mat[1, ] <- z0_vec
  z <- z0_vec
  
  # Expand the rate matrix
  R_prime <- matrix(rep(0, nevents * ntypes), c(nevents, ntypes))
  R_prime[cbind(1:nevents, p_vec)] <- r_vec
  
  while (t < tf)
  {
    event_rates <- R_prime %*% z  #multiply rates by pop sizes
    lambda <- sum(event_rates)  #combined rate at which stuff is happening
    if (lambda == 0)
    {
      # extinction
      zt_mat[times[times > t], ] <- rep(tail(zt_mat[rowSums(zt_mat) >= 0, ], 
                                         1), sum(times > t))  #fill the rest of the vector with last  datapoint
      return(out)
    }
    
    dt <- rexp(1, lambda)  #time until SOMETHING happends, don't know what yet
    
    event <- which.max(rmultinom(1, 1, event_rates))  #choose which event happened according to probabilities
    if (event > nevents)
    {
      print(event_rates)
      print(event)
    }
    increment <- e_mat[event, ]
    z <- z + increment  #add in the new offspring
    
    parent <- p_vec[event]
    z[parent] <- z[parent] - 1  #kill the parent
    
    times_to_update <- (t < times) & (t + dt >= times)
    zt_mat[times_to_update, ] <- rep(z, sum(times_to_update))
    t <- t + dt
  }
  
  return(zt_mat)
}

# simulating multiple branching processes where process rates vary as a
# function of dependent variables c_mat = dependent variable matrix.
# Dimensions c x q where q is nuber of dependent variables functions =
# vector of functions that calculate each model rate based on dependent
# vars. Dimensions m x 1 reps = vector of times to replicate each
# distinct condition. Dimensions c x 1

#' @export
bpsims <- function(model, theta, z0_vec, times, reps, c_mat = NA)
{
  if ((model$ndep > 0) && (is.na(c_mat) || ncol(c_mat) < model$ndep))
  {
    stop("Not enough dependent variables were provided for the model!")
  }
  if ((model$ndep == 0) && is.na(c_mat))
  {
    r_vec <- rep(0, ncol(model$e_mat))
    for (j in 1:nrow(model$e_mat))
    {
      r_vec[j] <- eval(model$func_deps[[j]], envir = list(c = theta))
    }
    x <- replicate(reps, bp(model$e_mat, r_vec, model$p_vec, z0_vec, times))
    z <- data.frame()
    for (i in 1:reps)
    {
      if (ncol(model$e_mat) == 1)
      {
        pop <- matrix(c(z0_vec, x[, , i]), ncol = 1)
      } else
      {
        pop <- matrix(rbind(z0_vec, x[, , i]), ncol = ncol(model$e_mat))
      }
      z <- rbind(z, data.frame(cbind(times = c(0, times), rep = i, 
                                     data.frame(pop))))
    }
    return(z)
    
  } else
  {
    if (length(reps) != nrow(c_mat))
    {
      stop("reps should be a vector with a different number of replications for each dependent variable condition")
    }
    z <- data.frame()
    r_vec <- 1
    for (i in 1:nrow(c_mat))
    {
      r_vec <- rep(0, ncol(model$e_mat))
      for (j in 1:nrow(model$e_mat))
      {
        r_vec[j] <- eval(model$func_deps[[j]], envir = list(c = theta, 
                                                        x = c_mat[i, ]))
      }
      x <- replicate(reps[i], bp(model$e_mat, r_vec, model$p_vec, z0_vec, times))
      for (k in 1:reps[i])
      {
        if (ncol(e_mat) == 1)
        {
          pop <- matrix(c(z0_vec, x[, , k]), ncol = 1)
        } else
        {
          pop <- matrix(rbind(z0_vec, x[, , k]), ncol = ncol(e_mat))
        }
        z <- rbind(z, data.frame(cbind(times = c(0, times), rep = r_vec, 
                                       variable_state = i, dep = c_mat[i, ], data.frame(pop))))
        r_vec <- r_vec + 1
      }
    }
  }
  
  return(z)
}


