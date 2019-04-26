library(dplyr)
library(ggplot2)

#Simulate multi-type Markov branching process
bp = function(B, D, Z0, times){#B = birth rate matrix, D =  death rates, Z0 = inital population vector, times = vector of timepoints to record
  tf = times[length(times)]
  d = length(D)
  Zt <- matrix(0, ncol = 2, nrow = length(times)) #store populations at different times
  t = 0
  Zt[1,] = Z0
  Z = Z0
  while(t < tf){
    eventRates = c(matrix(B*Z,d*d,1), D*Z) #rates at which all births and deaths are occuring given current population 
    lambda = sum(eventRates)
    if(lambda == 0){ #extinction
      out = list(pop = pop,time = times)
      return(out)
    }
    
    dt = rexp(1,lambda) #time until SOMETHING happends, don't know what yet

    event = which.max(rmultinom(1, 1, eventRates)) #choose which event happened according to probabilities
    
    if(event > d*d){ #event was a death
      Z[event - d*d] = Z[event - d*d]-1
    }
    else{ #event was a birth
      born = ceiling(event/d)
      Z[born] = Z[born]+1
    }
    
    times_to_update <- (t < times) & (t + dt >= times)
    Zt[times_to_update,] <- rep(Z, sum(times_to_update))
    t <- t + dt
  }
  
  return(Zt)
}


bpsims <- function(B, D, Z0, times, reps)
{
  x <- replicate(reps, bp(B, D, Z0, times))
  Z <- data.frame()
  for(i in 1:reps)
  {
    Z <- rbind(Z, data.frame(cbind(times, rep = i, data.frame(x[,,i]))))
  }
  return(Z)
}


library(expm)
B = matrix(c(.2,.02,.02,.2),ncol=2)# birth rate matrix
D = c(.15,.15) #death rate vector
d = length(D) # number of types
Z0 = c(1000,500) # initial population vector
Tf = 1 #final simulation timepoint

b = B/(rowSums(B) + D) #single individual offspring mean matrix
diag(b) = sapply(1:d, function(i){(B[i,i] + sum(B[i,]))/(sum(B[i,]) + D[i])})

#Diff EQ: M(t) = exp(At)
A =  B - diag(d)*(D)
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
dt = array(out[2,-1],c(d,d*d)) #second moment array

Mu = t(M)%*%Z0 #final mean population matrix
Sigma = matrix(Z0%*%dt,c(d,d*d)) #final population covariance matrix

for(i in 1:d){
  for(j in 1:d){
    Sigma[i,j] =  Sigma[i,j] - Z0%*%(M[,i]*M[,j])
  }
}

X = bpsims(B,D,Z0,times = times, reps = 50)
X <- X %>% group_by(rep)
names(X)[3:4] <- c("t1_cells", "t2_cells")
X <- X %>% group_by(rep) %>% mutate(t1_cells_prev = lag(t1_cells), t2_cells_prev = lag(t2_cells),
                                        t1_cells_0 = first(t1_cells), t2_cells_0 = first(t2_cells),
                                        dtimes = times - lag(times))

# Remove NA values
X <- X %>% filter(!is.na(t1_cells_prev))

pop_vec = cbind(X$t1_cells,X$t2_cells)
init_pop = cbind(X$t1_cells_prev,X$t2_cells_prev)
model_dtimes = sort(unique(times[-1] - times[-length(times)]))

mean(X[X$times == 1,]$t1_cells)
mean(X[X$times == 1,]$t2_cells)
var(X[X$times == 1,]$t1_cells)
var(X[X$times == 1,]$t2_cells)
cov(X[X$times == 1,]$t1_cells,X[X$times == 1,]$t2_cells)

library(rstan)
options(mc.cores = parallel::detectCores())
stan_dat <- list(d = d, N = nrow(pop_vec),  pop_vec = pop_vec, init_pop = init_pop, Tn = 1, times = as.array(model_dtimes), timesIdx = as.numeric(factor(X$dtimes)))
fit <- stan_model(file = "multitype_birth_death.stan")
fit.data <- sampling(fit, data = stan_dat, control = list(adapt_delta = 0.8))
s = extract(fit.data)


