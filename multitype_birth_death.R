library(dplyr)
library(ggplot2)

#Simulate multi-type Markov branching process
#B = birth rate matrix. Dimensions: N x (d+2) where N can be any number of offspring combinations with nonzero probability
#B[i, 1..d] contains the number of offspring of each type, B[i,d+1] contains ancestor type. B[i,d+2] is the rate at which 
                                                                                        #that combination will be born

bp = function(B, Z0, times){#D =  death rates, Z0 = inital population vector, times = vector of timepoints to record
  tf = times[length(times)]
  d = length(Z0)
  nEvents = nrow(B)
  Zt <- matrix(0, ncol = 2, nrow = length(times)) #store populations at different times
  t = 0
  Zt[1,] = Z0
  Z = Z0
  while(t < tf){
    eventRates = B[,d+2]*Z[B[,d+1]]
    lambda = sum(eventRates) #combined rate at which stuff is happening
    if(lambda == 0){ #extinction
      out = list(pop = pop,time = times)
      return(out)
    }
    
    dt = rexp(1,lambda) #time until SOMETHING happends, don't know what yet

    event = which.max(rmultinom(1, 1, eventRates)) #choose which event happened according to probabilities
    increment = B[event, 1:d]
    Z = Z + increment #add in the new offspring
    
    parent = B[event, d+1]
    Z[parent] = Z[parent] - 1 #kill the parent
    
    times_to_update <- (t < times) & (t + dt >= times)
    Zt[times_to_update,] <- rep(Z, sum(times_to_update))
    t <- t + dt
  }
  
  return(Zt)
}


bpsims <- function(B, Z0, times, reps)
{
  x <- replicate(reps, bp(B, Z0, times))
  Z <- data.frame()
  for(i in 1:reps)
  {
    Z <- rbind(Z, data.frame(cbind(times, rep = i, data.frame(x[,,i]))))
  }
  return(Z)
}


library(expm)
B2 = rbind(c(2, 0, 1, .2),c(1,1,1,.02),c(1,1,2,.02),c(0,2,2,.2), c(0,0,1,.15),c(0,0,2,.15))
B = matrix(c(.2,.02,.02,.2),ncol=2)# birth rate matrix
D = c(.15,.15)
d = length(D) # number of types
Z0 = c(100,50) # initial population vector
Tf = 5 #final simulation timepoint

b = matrix(rep(0,d*d),d,d)
lamb = rep(0,d)
for(i in 1:d){
  b[i,] = colSums(B2[B2[,d+1] == i, 1:d]*B2[B2[,d+1] == i,d+2])/sum(B2[B2[,d+1] == i,d+2]) #weighted average of columns gives expected progeny
  lamb[i] = sum(B2[B2[,d+1] == i,d+2])
}
#Diff EQ: M(t) = exp(At)
A =  lamb*(b - diag(d))
M = expm(A*Tf)

C = array(rep(0,d**3), c(d, d, d)); #matrix of second derivatives of offspring PGF
for(i in 1:d){
  C[i,,i] = b[i,]
  C[,i,i] = b[i,]
  diag(C[,,i]) = 0
  C[i,i,i] = (3*B[i,i] + sum(B[i,]))/(sum(B[i,]) + D[i]) - b[i,i] #- b[i,i]**2
  
}

library(deSolve)

second_moment_de = function(t, state, parameters){
  mt = expm(A*t)
  Beta = matrix(rep(0,d**3), nrow = d, ncol = d*d); #Beta array for second moment ODE
  for(i in 1:d){
    lamb = sum(B[i,]) + D[i]
    Beta[i,] = c(lamb*t(mt)%*%C[,,i]%*%mt) #vectorized computation of Beta
  }
  x = matrix(state,c(d,d*d))
  x_prime = c(A%*%x + Beta)
  return(list(x_prime))
}

init_state = c(1,0,0,0,0,0,0,1) #inital values of second moments
times = seq(0,5)

out = ode(y = init_state, times, func = second_moment_de, parms = 0)
dt = array(out[nrow(out),-1],c(d,d*d)) #second moment array

Mu = t(M)%*%Z0 #final mean population matrix
Sigma = matrix(Z0%*%dt,c(d,d*d)) #final population covariance matrix

for(i in 1:d){
  for(j in 1:d){
    Sigma[i,j] =  Sigma[i,j] - Z0%*%(M[,i]*M[,j])
  }
}

X = bpsims(B2, Z0,times = times, reps = 1000)
X <- X %>% group_by(rep)
names(X)[3:4] <- c("t1_cells", "t2_cells")
X <- X %>% group_by(rep) %>% mutate(t1_cells_prev = lag(t1_cells), t2_cells_prev = lag(t2_cells),
                                        t1_cells_0 = first(t1_cells), t2_cells_0 = first(t2_cells),
                                        dtimes = times - lag(times))

# Remove NA values
X <- X %>% filter(!is.na(t1_cells_prev))

#pop_vec = cbind(X$t1_cells,X$t2_cells)
#init_pop = cbind(X$t1_cells_prev,X$t2_cells_prev)
#model_dtimes = sort(unique(times[-1] - times[-length(times)]))

#mean(X[X$times == 1,]$t1_cells)
#mean(X[X$times == 1,]$t2_cells)
#var(X[X$times == 1,]$t1_cells)
#var(X[X$times == 1,]$t2_cells)
#cov(X[X$times == 1,]$t1_cells,X[X$times == 1,]$t2_cells)

#library(rstan)
#options(mc.cores = parallel::detectCores())
#stan_dat <- list(d = d, N = nrow(pop_vec),  pop_vec = pop_vec, init_pop = init_pop, Tn = 1, times = as.array(model_dtimes), timesIdx = as.numeric(factor(X$dtimes)))
#fit <- stan_model(file = "multitype_birth_death.stan")
#fit.data <- sampling(fit, data = stan_dat, control = list(adapt_delta = 0.8))
#s = extract(fit.data)


