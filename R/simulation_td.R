#' Simulate multi-type Markov branching process
#'
#' @param e_mat The matrix of birth events that can occur in the branching process. Dimensions \code{nevents} x \code{ntypes}
#' @param p_vec A vector containing the parent type for each of the birth events in \code{e_mat}. Dimensions \code{nevents} x 1
#' @param r_func A function of time comupting the vector containing the rate at which each of the birth events in \code{e_mat} occurs at time t.
#' @param r_func_ub A function of time computing upper bounds for the event rates after time t
#' @param z0_vec The initial population vector at time 0. Dimensions \code{ntypes} x 1
#' @param times The timepoints at which to record the population
#' @return A matrix with the population vectors at each timepoint. Dimensions \code{ntimes} x \code{ntypes}
bp_td <- function(e_mat, r_fun, r_fun_ub, p_vec, z0_vec, times)
{
  # z0_vec = inital population vector, times = vector of timepoints to record
  tf <- times[length(times)]
  ntypes <- length(z0_vec)
  nevents <- nrow(e_mat)
  zt_mat <- matrix(0, ncol = ntypes, nrow = length(times))  #store populations at different times
  t <- 0
  zt_mat[1, ] <- z0_vec
  z <- z0_vec
  
  while (t < tf & sum(z) > 0)
  {
    #Compute the time of the next event by simulating a NHPP arrival with the thinning method
    accept <- F
    t_next <- t
    while(!accept && t_next < tf){
      total_ub <- sum(r_fun_ub(t_next)*z[p_vec])
      t_next <- t_next + rexp(1,total_ub)
      y <- runif(1,0, sum(r_fun_ub(t_next)*z[p_vec]))
      if(y < sum(r_fun(t_next)*z[p_vec])){
        accept <- T
      }
    }
    
    event_rates <- r_fun(t_next)*z[p_vec]  #multiply rates by pop sizes
    
    event <- which.max(rmultinom(1, 1, event_rates))  #choose which event happened according to probabilities
    increment <- e_mat[event, ]
    z <- z + increment  #add in the new offspring
    
    parent <- p_vec[event]
    z[parent] <- z[parent] - 1  #kill the parent
    
    times_to_update <- (t < times) & (t_next >= times)
    zt_mat[times_to_update, ] <- rep(z, sum(times_to_update))
    t <- t_next
  }
  
  return(zt_mat)
}

#' Simulate multi-type Markov branching process where process rates vary as a function of dependent variables
#'
#' @param model An object of type \code{bpmodel} representing the model to be simulated
#' @param theta A vector of dependent variables to calculate the rates of each birth venet in \code{model}. Dimensions \code{ndep_levels} x \code{ndep}
#' @param z0_vec The initial population vector at time 0. Dimensions \code{ntypes} x 1
#' @param times The timepoints at which to record the population
#' @param reps a vector containing the number of times to replicate each dependent variables condidion. Dimensions \code{ndep_levels} x 1
#' 
#' @return A matrix with the population vectors at each timepoint. Dimensions \code{ntimes} x \code{ntypes}
#'
#' @export
bpsims_td <- function(model, theta, z0_vec, times, reps, c_mat = NA)
{
  if(attr(model, "class") != "bp_model"){
    stop("model must be a bp_model object!")
  }
  if ((model$ndep > 0) && (is.na(c_mat) || ncol(c_mat) != model$ndep))
  {
    stop("Incorrect number of dependent variables were provided for the model!")
  }
  if(length(theta) != model$nparams){
    stop("Wrong number of parameters were provided for the model!")
  }
  if(length(z0_vec) != ncol(model$e_mat)){
    stop("Initial popualation vector has wrong size!")
  }
  if ((model$ndep == 0) && is.na(c_mat))
  {
    r_vec <- rep(0, ncol(model$e_mat))
    for (j in 1:nrow(model$e_mat))
    {
      r_vec[j] <- eval(model$func_deps[[j]], envir = list(c = theta))
    }
    x <- array(replicate(reps, bp(model$e_mat, r_vec, model$p_vec, z0_vec, times)),c(length(times),length(z0_vec), reps))
    z <- data.frame()
    for (i in 1:reps)
    {
      pop <- matrix(rbind(x[, , i]), ncol = ncol(model$e_mat))
      z <- rbind(z, data.frame(cbind(times = times, rep = i, 
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
    curr_rep <- 1
    for (i in 1:nrow(c_mat))
    {
      r_vec <- rep(0, ncol(model$e_mat))
      for (j in 1:nrow(model$e_mat))
      {
        r_vec[j] <- eval(model$func_deps[[j]], envir = list(c = theta, 
                                                            x = c_mat[i, ]))
      }
      x <- array(replicate(reps[i], bp(model$e_mat, r_vec, model$p_vec, z0_vec, times)),c(length(times),length(z0_vec), reps[i]))
      for (k in 1:reps[i])
      {
        
        pop <- matrix(x[, , k], ncol = ncol(c_mat))
        dep = matrix(rep(c_mat[i, ], length(times)), ncol = ncol(c_mat))
        z <- rbind(z, data.frame(cbind(times = times, rep = curr_rep, 
                                       variable_state = i, data.frame(dep), data.frame(pop))))
        curr_rep <- curr_rep + 1
      }
    }
  }
  
  return(z)
}