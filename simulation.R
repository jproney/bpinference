library(dplyr)
library(ggplot2)

#Simulate multi-type Markov branching process
#E = event matrix. E[i] has number of each offspring produced in birth event i
#R = rate vector. Dimensions: N x 1. R[i] has rate of bith event i.
#P = parent vector. Dimensions: N x 1. P[i] has parent type of birth event i.

bp = function(E, R, P, Z0, times){#Z0 = inital population vector, times = vector of timepoints to record
  tf = times[length(times)]
  d = length(Z0)
  nEvents = nrow(E)
  Zt <- matrix(0, ncol = d, nrow = length(times)) #store populations at different times
  t = 0
  Zt[1,] = Z0
  Z = Z0
  
  #Expand the rate matrix
  R_prime = matrix(rep(0,nEvents*d), c(nEvents,d))
  R_prime[cbind(1:nEvents, P)] = R

  while(t < tf){
    eventRates = R_prime%*%Z #multiply rates by pop sizes
    lambda = sum(eventRates) #combined rate at which stuff is happening
    if(lambda == 0){ #extinction
      Zt[times[times > t],] = rep(tail(Zt[rowSums(Zt) >= 0,],1), sum(times > t)) #fill the rest of the vector with last  datapoint
      return(out)
    }
    
    dt = rexp(1,lambda) #time until SOMETHING happends, don't know what yet

    event = which.max(rmultinom(1, 1, eventRates)) #choose which event happened according to probabilities
    if(event > nEvents){
      print(eventRates)
      print(event)
    }
    increment = E[event,]
    Z = Z + increment #add in the new offspring
    
    parent = P[event]
    Z[parent] = Z[parent] - 1 #kill the parent
    
    times_to_update <- (t < times) & (t + dt >= times)
    Zt[times_to_update,] <- rep(Z, sum(times_to_update))
    t <- t + dt
  }
  
  return(Zt)
}

# simulating multiple branching processes where process rates vary as a function of dependent variables
# C = dependent variable matrix. Dimensions c x q where q is nuber of dependent variables
# functions = vector of functions that calculate each model rate based on dependent vars. Dimensions m x 1
# reps = vector of times to replicate each distinct condition. Dimensions c x 1
bpsims = function(model, theta, Z0, times, reps, C = NA){
  if((model$nDep > 0) && (is.na(C) || ncol(C) < model$nDep)){
    stop("Not enough dependent variables were provided for the model!")
  }
  if((model$nDep == 0) && is.na(C)){
    R = rep(0,ncol(model$E))
    for(j in 1:nrow(model$E)){
      R[j] = eval(model$func_deps[[j]], envir = list(c = theta))
    }
    x <- replicate(reps, bp(model$E, R, model$P, Z0, times)) 
    Z <- data.frame()
    times = c(0,times)
    for(i in 1:reps)
    {
      pop = matrix(rbind(Z0,x[,,i]),ncol=ncol(E))
      Z <- rbind(Z, data.frame(cbind(times, rep = i, data.frame(pop))))
    }
    return(Z)
    
  }
  else{
    if(length(reps) != nrow(C)){
      stop("reps should be a vector with a different number of replications for each dependent variable condition")
    }
    Z <- data.frame()
    r = 1
    for(i in 1:nrow(C)){
      R = rep(0,ncol(model$E))
      for(j in 1:nrow(model$E)){
        R[j] = eval(model$func_deps[[j]], envir = list(c = theta, x = C[j,]))
      }
      x <- replicate(reps[i], bp(model$E, R, model$P, Z0, times)) 
      for(k in 1:reps[i])
      {
        Z <- rbind(Z, data.frame(cbind(rbind(0,times), rep = r, variable_state = i, rbind(Z0,data.frame(x[,,i])), C[i,])))
        r <- r + 1
      }
    }
  }
  
  return(Z)
}


