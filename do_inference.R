#wrapper script to aid in branching process inference

# pop_vec = population vectors at end of each run. Dimensions n x d
# init_pop = population at beginning of each run. DImensions n x d
# times = the length of time elapsing between each initial population and each final population. dimensions n x 1
# E = birth events matrix. Dimensions m x d
# P = parents matrix. Dimensions m x 1
# C = matrix of dependent variables that vary from run to run. Dimensions n x q
# functions = vector of function expressions specifying how each rate depends on the dependent variables in the C matrix
# Priors = priors for all of the parameters. Should be in the form of a list of vectors
# Priors should be a z-dimensional list, and each list entry should have a prior for that parameter
# if C or functions is left as NA, inference is performed directly on the rates as parameters

do_inference = function(final_pop, init_pop, E, P, times, priors, C = NA, functions = NA){
  m = nrow(E) #number of events
  d = ncol(E) #number of types
  n = nrow(final_pop) #number of datapoints
  times_unique = unique(times) #distinct timepoints
  l = length(times_unique) #number of distinct durations
  times_idx = match(times,times_unique) #index of time duration for each datapoint
  z = length(priors) #total number of parameters
  
  if(!is.na(C) && !is.na(functions)){
    C_unique = C[!duplicated(C),] #unique combinations of dependent variables
    c = nrow(C_unique) #number of distinct combinations of dependent variables
    q = ncol(C_unique) #number of different dependent variables
    var_idx = apply(C,1,function(r){which.min(abs(rowSums(sweep(C_unique,2,r))))})  #indices of unique dependent variable combinations
  }
  else{
    C_unique = matrix(0, 1,1)
    c = 1
    q = 1
    var_idx = rep(1, n)
    functions = rep(0,m)
    for(i in 1:m){
      functions[i] = parse(text=sprintf("c%d", i))
    }
  }
  
  library(rstan)
  options(mc.cores = parallel::detectCores())
  generate(functions, priors, "multitype_birth_death.stan") #generate the stan file
  stan_dat = list(d = d, m = m, n = n, l=l, c = c, q=q, z=z, E = E, P = P, 
                   pop_vec = final_pop, init_pop = init_pop,
                   times = array(times_unique,1), times_idx =  times_idx, 
                   function_var = C_unique, var_idx = var_idx)
  fit = stan_model(file = "multitype_birth_death.stan")
  fit.data = sampling(fit, data = stan_dat, control = list(adapt_delta = 0.8), chains = 4, refresh = 1, init="0")
  return(fit.data)
}