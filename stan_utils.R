#wrapper script to aid in branching process inference

# model = the branching process model being estimated
# pop_vec = population vectors at end of each run. Dimensions n x d
# init_pop = population at beginning of each run. DImensions n x d
# times = the length of time elapsing between each initial population and each final population. dimensions n x 1
# C = matrix of dependent variables that vary from run to run. Dimensions n x q
# functions = vector of function expressions specifying how each rate depends on the dependent variables in the C matrix
# Priors = priors for all of the parameters. Should be in the form of a list of vectors
# Priors should be a z-dimensional list, and each list entry should have a prior for that parameter
# if C or functions is left as NA, inference is performed directly on the rates as parameters

create_stan_data = function(model, final_pop, init_pop, times, priors, C = NA){
  m = nrow(model$E) #number of events
  d = ncol(model$E) #number of types
  n = nrow(final_pop) #number of datapoints
  times_unique = unique(times) #distinct timepoints
  l = length(times_unique) #number of distinct durations
  times_idx = match(times,times_unique) #index of time duration for each datapoint
  z = model$nParams #total number of parameters
  if(model$nParams != length(priors)){
    stop("Incorrect number of priors for model!")
  }
  
  if(model$nDep > 0 && (is.na(C) || ncol(C) < mod$nDep)){
    stop("C does not contain enough dependent variables for the model")
  }
  if(model$nDep == 0 && is.na(C)){
    C_unique = matrix(0, 1,1)
    c = 1
    q = 1
    var_idx = rep(1, n)
  }
  else{
    C_unique = matrix(C[!duplicated(C),],ncol = ncol(C)) #unique combinations of dependent variables
    c = nrow(C_unique) #number of distinct combinations of dependent variables
    q = ncol(C_unique) #number of different dependent variables
    var_idx = apply(C,1,function(r){which.min(abs(rowSums(sweep(C_unique,2,r))))})  #indices of unique dependent variable combinations
  }

  library(rstan)
  generate(model, priors, "multitype_birth_death.stan") #generate the stan file
  stan_dat = list(d = d, m = m, n = n, l=l, c = c, q=q, z=z, E = E, P = P, 
                   pop_vec = final_pop, init_pop = init_pop,
                   times = array(times_unique,1), times_idx =  times_idx, 
                   function_var = C_unique, var_idx = var_idx)
  fit = stan_model(file = "multitype_birth_death.stan")
  return(list(model = fit, data = stan_dat))
}

# create Stan initialization list by selecting initial values uniformly at ranom.
# ranges is a z x 2 matrix, where z is the number of parameters. Each row has contains a lower and upper bound for initialization.
uniform_initialize = function(ranges, nchains){
  return(replicate(nchains, list(list(Theta = apply(ranges,1,function(s){runif(1,s[1],s[2])})))))
}

# Helper function for piping simulation output into inference engine
stan_data_from_simulation = function(X, model){
  d = ncol(model$E)
  for(i in 1:d){
    cellname = sprintf("t%d_cells", i)
    names(X)[2 + i] = cellname
    prevname = paste(cellname,"prev",sep="_")
    X = X %>% mutate(prev = lag(X[,cellname]))
    names(X)[2 + d + i] = prevname
  }
  return(X)
}
