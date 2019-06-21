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

create_stan_data = function(model, final_pop, init_pop, times, C = NA){
  m = nrow(model$E) #number of events
  d = ncol(model$E) #number of types
  n = nrow(final_pop) #number of datapoints
  times_unique = unique(times) #distinct timepoints
  l = length(times_unique) #number of distinct durations
  times_idx = match(times,times_unique) #index of time duration for each datapoint
  z = model$nParams #total number of parameters
  
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


  stan_dat = list(d = d, m = m, n = n, l=l, c = c, q=q, z=z, E = E, P = P, 
                   pop_vec = final_pop, init_pop = init_pop,
                   times = array(times_unique,1), times_idx =  times_idx, 
                   function_var = C_unique, var_idx = var_idx)
  return(stan_dat)
}

# create Stan initialization list by selecting initial values uniformly at ranom.
# ranges is a z x 2 matrix, where z is the number of parameters. Each row has contains a lower and upper bound for initialization.
uniform_initialize = function(ranges, nchains){
  return(replicate(nchains, list(list(Theta = apply(ranges,1,function(s){runif(1,s[1],s[2])})))))
}

# Helper function for piping simulation output into inference engine
stan_data_from_simulation = function(X, model){
  library(dplyr)
  d = ncol(model$E)
  offset = model$nDep + (model$nDep > 0)
  for(i in 1:d){
    cellname = sprintf("t%d_cells", i)
    names(X)[2 + offset + i] = cellname
    prevname = paste(cellname,"prev",sep="_")
    X = X %>% mutate(prev = lag(X[,cellname]))
    names(X)[2 + offset + d + i] = prevname
  }
  X = X %>% mutate(dtimes = times-lag(times))
  X <- X %>% filter(times!=0)
  init_pop = as.matrix(X[,(3+ offset + d):(2+ offset + 2*d)])
  final_pop = as.matrix(X[,(3 + offset):(2+ offset + d)])
  times = X$dtimes
  if(model$nDep > 0){
    return(create_stan_data(model, final_pop, init_pop, times, matrix(X[,4:(3+model$nDep)], ncol = model$nDep)))
  }
  return(create_stan_data(model, final_pop, init_pop, times))
}

#functional_deps is an array of strings encoding functions of  variable 'x1, x2, ...' and constant parameters 'c1','c2',...
#priors is list of with prior objects for all parameters in order of the rate they contribute to
#example prior object: p = list(name="gamma",params=c(1,1), init = c(0,3), bounds=c(0,3))
#filename is name of generated Stan file
generate = function(model, priors, filename){
  
  template =  readChar("stan_template.txt", file.info("stan_template.txt")$size) #load the template file
  nFunc = length(model$func_deps)
  funcs = rep(0,nFunc)
  
  dists = readLines("allowed_distributions.txt") #load distributions file
  nPri = length(priors)
  declStrs = rep(0,nPri)
  priorStrs = rep(0, nPri)
  
  for(i in 1:nFunc){
    # convert the expression
    exprn = model$func_deps[[i]]
    stanstr = exp_to_stan(exprn, model$nParams, mod$nDep)
    funcs[i] = sprintf("\t\tR[%d, P[%d]] = %s;\n",i,i,stanstr) 
  }
  
  for(p in 1:nPri){
    ind = which(!is.na(str_extract(dists, paste("^",priors[[p]]$name, " ",sep=""))))
    if(length(ind) > 0){
      nparm = strtoi(str_extract(dists[ind], "[1-9]+"))
      if(nparm != length(priors[[p]]$params)){
        stop("Incorrect number of parameters for prior distribution")
      }
      constraints = str_extract(dists[ind], "[_+]+")
      for(i in 1:nparm){
        if(priors[[p]]$params[i] <= 0 && substr(constraints,i,i) == "+"){
          stop("Negative parameter not allowed here!")
        }
      }
    }
    else{
      stop("Invalid prior name!")
    }
    priorStrs[p] = sprintf("\tTheta%d ~ %s(%s);\n", p, priors[[p]]$name, paste(priors[[p]]$params, collapse=","))
    if(!any(is.na(priors[[p]]$bounds))){
      declStrs[p] = sprintf("\treal<lower=%d, upper=%d> Theta%d;\n", priors[[p]]$bounds[1], priors[[p]]$bounds[2], p)
    }
  }
  write(sprintf(template, paste(declStrs, collapse=""), paste(funcs, collapse=""),  paste(priorStrs,collapse="")), filename)
}

# takes simple mathematical R expressions and turns them into a specific type of Stan code to fill in the template
exp_to_stan = function(exprn, maxParams, maxDep){
  if(is.atomic(exprn) || is.name(exprn)){
    first = exprn
  }
  else{
    first = exprn[[1]]
  }
  
  if(is.atomic(first)){
    return(first)
  }
  if(deparse(first) %in% c('+','-','/','*','^','exp', 'log')){ #allowed operations
    if(length(exprn) > 2){
      return(paste("(",exp_to_stan(exprn[[2]], maxParams, maxDep),")", deparse(first), "(", exp_to_stan(exprn[[3]], maxParams, maxDep), ")", sep=" "))
    }
    else{
      return(paste(deparse(first),'(', exp_to_stan(exprn[[2]], maxParams, maxDep), ')', sep=""))
    }
  }
  if(deparse(first) == '('){
    return(exp_to_stan(exprn[[2]], maxParams, maxDep))
  }
  if(deparse(first) == '['){
    if(!is.na(str_extract(deparse(exprn),"^c\\[[1-9]+\\]$"))){ # function parameters
      s = str_extract(deparse(exprn),"^c\\[[1-9]+\\]$")
      num = strtoi(substr(s,3,nchar(s)-1))
      if(num > maxParams){
        stop(paste("Parameter",s,"goes beyond the number of paramters specified.",sep=" "))
      }
      return(sprintf("Theta%d", num))
    }
    if(!is.na(str_extract(deparse(exprn),"^x\\[[1-9]+\\]$"))){ # variables
      s = str_extract(deparse(exprn),"^x\\[[1-9]+\\]$")
      num = strtoi(substr(s,3,nchar(s)-1))
      if(num > maxDep){
        stop(paste("Variable",s,"goes beyond the number of dependent variables specified.",sep=" "))
      }
      return(sprintf("function_var[i,%d]", num))
    }
  }
  else{
    stop("Invalid expression!")
  }
}

check_valid = function(exprn, maxParams, maxDep){
  if(is.atomic(exprn) || is.name(exprn)){
    first = exprn
  }
  else{
    first = exprn[[1]]
  }
  
  if(is.atomic(first)){
    return()
  }
  if(deparse(first) %in% c('+','-','/','*','^','exp', 'log')){ #allowed operations
    if(length(exprn) > 2){
      check_valid(exprn[[2]], maxParams, maxDep)
      check_valid(exprn[[3]], maxParams, maxDep)
      return()
    }
    else{
      check_valid(exprn[[2]], maxParams, maxDep)
      return();
    }
  }
  if(deparse(first) == '('){
    check_valid(exprn[[2]], maxParams, maxDep)
    return()
  }
  if(deparse(first) == '['){
    if(!is.na(str_extract(deparse(exprn),"^c\\[[1-9]+\\]$"))){ # function parameters
      s = str_extract(deparse(exprn),"^c\\[[1-9]+\\]$")
      num = strtoi(substr(s,3,nchar(s)-1))
      if(num > maxParams){
        stop(paste("Parameter",s,"goes beyond the number of paramters specified.",sep=" "))
      }
      return()
    }
    if(!is.na(str_extract(deparse(exprn),"^x\\[[1-9]+\\]$"))){ # variables
      s = str_extract(deparse(exprn),"^x\\[[1-9]+\\]$")
      num = strtoi(substr(s,3,nchar(s)-1))
      if(num > maxDep){
        stop(paste("Variable",s,"goes beyond the number of dependent variables specified.",sep=" "))
      }
      return()
    }
    stop(paste("Invalid expression:",deparse(exprn),sep=" "))
  }
  else{
    stop(paste("Invalid expression:",deparse(exprn),sep=" "))
  }
}
