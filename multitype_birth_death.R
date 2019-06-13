library(dplyr)
library(ggplot2)

#Simulate multi-type Markov branching process
#B = offspring matrix. Dimensions: N x d where N can be any number of offspring combinations with nonzero probability.
#E[i] has number of each offspring produced in birth event i
#R = rate matrix. Dimensions: N x 1. R[i] has rate of bith event i.
#P = parent matrix. Dimensions: N x 1. P[i] has parent type of birth event i.

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
      out = list(pop = pop,time = times)
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


bpsims <- function(E, R, P, Z0, times, reps)
{
  x <- replicate(reps, bp(E, R, P, Z0, times))
  Z <- data.frame()
  for(i in 1:reps)
  {
    Z <- rbind(Z, data.frame(cbind(times, rep = i, data.frame(x[,,i]))))
  }
  return(Z)
}




