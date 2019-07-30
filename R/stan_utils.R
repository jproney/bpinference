#' function to generate a stan input list from population data
#' 
#' @param model The branching process model being estimated
#' @param final_pop Population vectors at end of each run. Dimensions ndatapts x ntypes 
#' @param init_pop = population at beginning of each run.  Dimensions ndatapts x ntypes 
#' @param times = the length of time elapsing between each initial population and each final population. dimensions ndatapts x 1 
#' @param c_mat = matrix of dependent variables that vary from run to run. Dimensions ndatapts x ndep 
#' @param simple_bd = logical value indicating whether or not the model is a simple birth-death process. If TRUE, a simpler data list will be returned
#'
#' @return A data list to pass to the Stan sampling function 
#' @export
create_stan_data <- function(model, final_pop, init_pop, times, c_mat = NA, simple_bd = FALSE)
{
  if (is.vector(final_pop))
  {
    warning("final_pop is a vector, not a matrix. Converting to a one-column matrix")
    final_pop <- matrix(final_pop, ncol = 1)
  }
  if (is.vector(init_pop))
  {
    warning("init_pop is a vector, not a matrix. Converting to a one-column matrix")
    init_pop <- matrix(init_pop, ncol = 1)
  }
  if(nrow(final_pop) != nrow(init_pop) || nrow(init_pop) != length(times)){
    stop("times, initial populations, and final populations must all have same number of rows!")
  }
  nevents <- nrow(model$e_mat)  #number of events
  ntypes <- ncol(model$e_mat)  #number of types
  ndatapts <- nrow(final_pop)  #number of datapoints
  times_unique <- unique(times)  #distinct timepoints
  ntimes_unique <- length(times_unique)  #number of distinct durations
  times_idx <- match(times, times_unique)  #index of time duration for each datapoint
  nparams <- model$nparams  #total number of parameters
  
  if (model$ndep > 0 && (is.na(c_mat) || ncol(c_mat) < model$ndep))
  {
    stop("c_mat does not contain enough dependent variables for the model")
  }
  if (model$ndep == 0 && is.na(c_mat))
  { #fill c_mat with dummy data
    c_mat_unique <- matrix(0, 1, 1)
    ndep_levels <- 1
    ndep <- 1
    var_idx <- rep(1, ndatapts)
  } else
  {
    c_mat_unique <- matrix(c_mat[!duplicated(c_mat), ], ncol = ncol(c_mat))  #unique combinations of dependent variables
    ndep_levels <- nrow(c_mat_unique)  #number of distinct combinations of dependent variables
    ndep <- ncol(c_mat_unique)  #number of different dependent variables
    var_idx <- apply(c_mat, 1, function(r)
    {
      which.min(rowSums(abs(sweep(c_mat_unique, 2, r))))
    })  #indices of unique dependent variable combinations
  }
  
  if(!simple_bd){
    stan_dat <- list(ntypes = ntypes, nevents = nevents, ndatapts = ndatapts, ntimes_unique = ntimes_unique, ndep_levels = ndep_levels, ndep = ndep, nparams = nparams, e_mat = model$e_mat, 
      p_vec = model$p_vec, pop_vec = final_pop, init_pop = init_pop, times = array(times_unique,1), times_idx = times_idx, function_var = c_mat_unique, var_idx = var_idx)
  }
  else{
    stan_dat <- list(pop_vec = c(final_pop), init_pop = c(init_pop), times = times, function_var = c_mat_unique, var_idx = var_idx, 
                     ndep_levels = ndep_levels, ndep = ndep, nparams = nparams,ndatapts = ndatapts)
  }
  
  return(stan_dat)
}

#' create Stan initialization list by selecting initial values uniformly at ranom.  
#' @param range nparams x 2 matrix, where nparams is the number ofparameters. Each row has a lower and upper bound for initialization
#' @param nchains number of Markov chains for which to sample initial values
#' 
#' @return a stan initialization list
#' @export
uniform_initialize <- function(ranges, nchains)
{
  init_list <- list()
  names_list <- sapply(1:nrow(ranges), function(i){sprintf("Theta%d", i)})
  for(i in 1:nchains){
    init_list[[i]] <- as.list(apply(ranges, 1, function(s)
    {
      runif(1, s[1], s[2])
    }))
    names(init_list[[i]]) <- names_list 
    print(init_list)
  }
  return(init_list)
}

#' Helper function for piping simulation output into inference engine
#' @param sim_data simulation data as returned from \code{bpsims}
#' @param model the \code{bpmodel} object the produced the data
#' @param simple_bd = logical value indicating whether or not the model is a simple birth-death process. If TRUE, a simpler data list will be returned
#' 
#' @return A data list to pass to the Stan sampling function 
#' @export
stan_data_from_simulation <- function(sim_data, model, simple_bd = FALSE)
{
  ntypes <- ncol(model$e_mat)
  offset <- model$ndep + (model$ndep > 0)
  for (i in 1:ntypes)
  {
    cellname <- sprintf("t%d_cells", i)
    names(sim_data)[2 + offset + i] <- cellname
    prevname <- paste(cellname, "prev", sep = "_")
    sim_data <- dplyr::group_by(sim_data, rep)
    sim_data <- dplyr::mutate(sim_data, prev = dplyr::lag(!!sym(cellname)))
    names(sim_data)[2 + offset + ntypes + i] <- prevname
  }
  sim_data <- dplyr::mutate(sim_data,  dtimes = times - dplyr::lag(times))
  sim_data <- dplyr::filter(sim_data,  !is.na(t1_cells_prev))
  init_pop <- as.matrix(sim_data[, (3 + offset + ntypes):(2 + offset + 2 * ntypes)])
  final_pop <- as.matrix(sim_data[, (3 + offset):(2 + offset + ntypes)])
  times <- sim_data$dtimes
  if (model$ndep > 0)
  {
    return(create_stan_data(model, final_pop, init_pop, times, as.matrix(sim_data[, 
                                                                        4:(3 + model$ndep)], ncol = model$ndep), simple_bd = simple_bd))
  }
  return(create_stan_data(model, final_pop, init_pop, times, simple_bd = simple_bd))
}

#' Generate a Stan file for inferring the model parameters
#' @param model The \code{bpmodel} object from which to generate the Stan file
#' @param priors A list of prior distributions for each parameters in the mode. 
#' Each entry in the list should be a names list of the form list(name="normal",params=c(0,1),bounds=(-100,100))
#' @param filename The name of the file in which to store the model. Should end in ".stan"
#' 
#' @return a string containing the model code
#' @export
generate <- function(model, priors, filename = NA, simple_bd = FALSE, noise_model = F)
{
  if(is.null(attr(model, "class")) ||  attr(model, "class") != "bp_model"){
    stop("model must be a bp_model object!")
  }
  if(!is.na(filename) && (!is.vector(filename) || !is.character(filename))){
    stop("filename should be a character vector!")
  }
  nfunc <- length(model$func_deps)
  funcs <- rep(0, nfunc)
  
  nprior <- length(priors)
  declStrs <- rep(0, nprior)
  priorStrs <- rep(0, nprior)
  
  if(nprior != model$nparams){
    stop("incorrect number of priors for model!")
  }
  
  for (i in 1:nfunc)
  {
    # convert the expression
    exprn <- model$func_deps[[i]]
    stanstr <- exp_to_stan(exprn, model$nparams, model$ndep)
    if(!simple_bd){
      funcs[i] <- sprintf("\t\tr_mat[%d, p_vec[%d]] = %s;\n", i, i, stanstr)
    }else{
      funcs[i] <- sprintf("\t\tr_mat[%d, i] = %s;\n", i, stanstr)
    }
  }

  for (p in 1:nprior)
  {
    parse <- parse_prior(priors[[p]],p)
    priorStrs[p] <- parse[1]
    declStrs[p] <- parse[2]
  }
  if(simple_bd){
    if(!all(model$e_mat == matrix(c(2,0),ncol=1))){
      stop("model is not a simple birth-death process!")
    }
    model_str = sprintf(STAN_TEMPLATE_SIMPLE_BD, paste(declStrs, collapse = ""), paste(funcs, 
                                                                collapse = ""), paste(priorStrs, collapse = ""))
  }
  else{
    if(noise_model){
      if(ncol(model$e_mat) == 1){
        stop("observational noise models are not available for one-type processes")
      }
      model_str = sprintf(STAN_TEMPLATE_MULTITYPE_NOISE, paste(declStrs, collapse = ""), paste(funcs, 
                                                                collapse = ""), paste(priorStrs, collapse = ""))
    }
    else{
      model_str = sprintf(STAN_TEMPLATE_MULTITYPE, paste(declStrs, collapse = ""), paste(funcs, 
                                                                collapse = ""), paste(priorStrs, collapse = ""))
    }
  }
  if(!is.na(filename)){
    write(model_str, filename)
  }
  return(model_str)
}

#' takes simple mathematical R expressions and turns them into a specific type of Stan code to fill in the template
#' @param exprn the expression to parse
#' @param max_params the number of parameters in the model, which determines the maximum parameter index that can be referenced in the expression
#' @param max_dep the number of dependent variables in the model, which determines the amximum dependent variables index that can be referenced in the model
#' 
#' @return A string with the Stan code derived from the R expressions
exp_to_stan <- function(exprn, max_params, max_dep)
{
  if (is.atomic(exprn) || is.name(exprn))
  {
    first <- exprn
  } else
  {
    first <- exprn[[1]]
  }
  
  if (is.atomic(first) && is.numeric(first) && !is.array(first))
  {
    return(first)
  }
  if (deparse(first) %in% c("+", "-", "/", "*", "^"))
  {
    # allowed operations
    if (length(exprn) > 2)
    {
      return(paste("(", exp_to_stan(exprn[[2]], max_params, max_dep), 
                   ")", deparse(first), "(", exp_to_stan(exprn[[3]], max_params, 
                                                         max_dep), ")", sep = " "))
    } else
    {
      return(paste(deparse(first), "(", exp_to_stan(exprn[[2]], max_params, 
                                                    max_dep), ")", sep = ""))
    }
  }
  if(deparse(first) %in% c("exp", "log"))
  {
    if (length(exprn) > 2)
    {
      stop("log and exp take only one parameter")  
    } else
    {
      return(paste(deparse(first), "(", exp_to_stan(exprn[[2]], max_params, 
                                                    max_dep), ")", sep = ""))
    }
      
  }
  if (deparse(first) == "(")
  {
    return(exp_to_stan(exprn[[2]], max_params, max_dep))
  }
  if (deparse(first) == "[")
  {
    if (!is.na(stringr::str_extract(deparse(exprn), "^c\\[[0-9]+\\]$")))
    {
      # function parameters
      s <- stringr::str_extract(deparse(exprn), "^c\\[[0-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_params || num <= 0)
      {
        stop(paste("Parameter", s, "goes beyond the number of parameters specified.", 
                   sep = " "))
      }
      return(sprintf("Theta%d", num))
    }
    if (!is.na(stringr::str_extract(deparse(exprn), "^x\\[[0-9]+\\]$")))
    {
      # variables
      s <- stringr::str_extract(deparse(exprn), "^x\\[[0-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_dep || num <= 0)
      {
        stop(paste("Variable", s, "goes beyond the number of dependent variables specified.", 
                   sep = " "))
      }
      return(sprintf("function_var[i,%d]", num))
    }
    stop(paste("Invalid expression:", deparse(exprn), sep = " "))
  } 
  else
  {
    stop(paste("Invalid expression:", deparse(exprn), sep = " "))
  }
}

#' takes simple mathematical R expressions and determines whether they can be parsed into stan code
#' @param exprn the expression to parse
#' @param max_params the number of parameters in the model, which determines the maximum parameter index that can be referenced in the expression
#' @param max_dep the number of dependent variables in the model, which determines the amximum dependent variables index that can be referenced in the model
#' 
#' @return Nothing if the expression is valid. Otherwise an exception is thrown.
check_valid <- function(exprn, max_params, max_dep)
{
  if (is.atomic(exprn) || is.name(exprn))
  {
    first <- exprn
  } else
  {
    first <- exprn[[1]]
  }
  
  if (is.atomic(first) && is.numeric(first) && !is.array(first))
  {
    return()
  }
  if (deparse(first) %in% c("+", "-", "/", "*", "^"))
  {
    # allowed operations
    if (length(exprn) > 2)
    {
      check_valid(exprn[[2]], max_params, max_dep)
      check_valid(exprn[[3]], max_params, max_dep)
      return()
    } else
    {
      check_valid(exprn[[2]], max_params, max_dep)
      return()
    }
  }
  if(deparse(first) %in% c("exp", "log"))
  {
    if (length(exprn) > 2)
    {
      stop("log and exp take only one parameter")  
    } else
    {
      check_valid(exprn[[2]], max_params, max_dep)
      return()
    }
  }
  if (deparse(first) == "(")
  {
    check_valid(exprn[[2]], max_params, max_dep)
    return()
  }
  if (deparse(first) == "[")
  {
    if (!is.na(stringr::str_extract(deparse(exprn), "^c\\[[0-9]+\\]$")))
    {
      # function parameters
      s <- stringr::str_extract(deparse(exprn), "^c\\[[0-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_params || num <= 0)
      {
        stop(paste("Parameter", s, "goes beyond the number of parameters specified.", 
                   sep = " "))
      }
      return()
    }
    if (!is.na(stringr::str_extract(deparse(exprn), "^x\\[[0-9]+\\]$")))
    {
      # variables
      s <- stringr::str_extract(deparse(exprn), "^x\\[[0-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_dep || num <= 0)
      {
        stop(paste("Variable", s, "goes beyond the number of dependent variables specified.", 
                   sep = " "))
      }
      return()
    }
    stop(paste("Invalid expression:", deparse(exprn), sep = " "))
  } else
  {
    stop(paste("Invalid expression:", deparse(exprn), sep = " "))
  }
}

#' parses prior object into a prior statement that can be inserted into a Stan model (e.g. Theta4 ~ normal(0,1))
#' @param prior the \code{prior_dist} object to be parsed
#' @param k the number parameter corresponding to the given prior.
#' 
#' @return 2 strings: A parameter declaration the goes in the \code{parameters} block of the Stan model, 
#' and a prior probability statement that goes in the \code{model} block of the Stan model
parse_prior = function(prior, k){
  ind <- which(!is.na(stringr::str_extract(ALLOWED_DISTRIBUTIONS, paste("^", prior$name, 
                                                                          " ", sep = ""))))
  if (length(ind) > 0)
  {
    nparm <- strtoi(stringr::str_extract(ALLOWED_DISTRIBUTIONS[ind], "[1-9]+"))
    if (nparm != length(prior$params))
    {
      stop("Incorrect number of parameters for prior distribution")
    }
    constraints <- stringr::str_extract(ALLOWED_DISTRIBUTIONS[ind], "[_+]+")
    for (i in 1:nparm)
    {
      if (prior$params[i] <= 0 && substr(constraints, i, 
                                                 i) == "+")
      {
        stop("Positive parameter required here!")
      }
    }
  } else
  {
    stop("Invalid prior name!")
  }
  priorStr <- sprintf("\tTheta%d ~ %s(%s);\n", k, prior$name, 
                          paste(prior$params, collapse = ","))
  if (!any(is.na(prior$bounds)))
  {
    if(!is.numeric(prior$bounds) || length(prior$bounds) != 2){
      stop("bounds should be a list of 2 numbers!")
    }
    if(prior$bounds[1] >= prior$bounds[2]){
      stop("Error: lower bounder greater or equal to upper bound!")
    }
    declStr <- sprintf("\treal<lower=%s, upper=%s> Theta%d;\n", 
                           paste(prior$bounds[1]), paste(prior$bounds[2]), k)
  }
  else{
    declStr <- sprintf("\treal Theta%d;\n", k) 
  }
  return(c(priorStr, declStr))
}

#' function to make sure a \code{prior_dist} object actually encodes a vaild prior distribution
#' 
#' @param prior, the \code{prior_dist} object being validated
#' @return the \code{prior_dist} if the prior is valid. Will throw an error otherwise.
validate_prior_dist <- function(prior){
  if(attr(prior,'class') != "prior_dist"){
    stop("validate_prior_dist required an object of type prior_dist")
  }
  data <- unclass(prior)
  ind <- which(!is.na(stringr::str_extract(ALLOWED_DISTRIBUTIONS, paste("^", data$name, 
                                                                        " ", sep = ""))))
  if (length(ind) > 0)
  {
    nparm <- strtoi(stringr::str_extract(ALLOWED_DISTRIBUTIONS[ind], "[1-9]+"))
    if (nparm != length(data$params))
    {
      
      stop("Incorrect number of parameters for prior distribution")
    }
    constraints <- stringr::str_extract(ALLOWED_DISTRIBUTIONS[ind], "[_+]+")
    for (i in 1:nparm)
    {
      if (data$params[i] <= 0 && substr(constraints, i, 
                                         i) == "+")
      {
        stop("Positive parameter required here!")
      }
    }
    if(data$name == "uniform" && data$params[2] <= data$params[1]){
      stop("upper bound must be greater than lower in uniform distribution!")
    }
  } else
  {
    stop("Invalid prior name!")
  }
  if (!any(is.na(data$bounds)))
  {
    if(!is.numeric(data$bounds) || length(data$bounds) != 2){
      stop("bounds should be a list of 2 numbers!")
    }
    if(data$bounds[1] >= data$bounds[2]){
      stop("Error: lower bounder greater or equal to upper bound!")
    }
  }
  return(prior)
}

#' Creates a new \code{prior_dist} object
#' 
#' @param name a string with the name of the distribution (e.g. 'normal'). Must match a distribution in \code{ALLOWED_DISTRIBUTIONS} list
#' @param params parameters for the prior distribution. Should match the number and constrains imposed by the distribution name.
#' @param bounds a lower and upper bound for the parameter. This restricts the range in which sampling occurs. If left as NA, no bounds are imposed
#' @return a new \code{prior_dist} object
new_prior_dist <- function(name, params, bounds = NA){
  return(structure(list(name = name, params = params, bounds = bounds), class = "prior_dist"))
}

#' User-facing constructor for new \code{prior_dist} objects
#' 
#' @param name a string with the name of the distribution (e.g. 'normal'). Must match a distribution in \code{ALLOWED_DISTRIBUTIONS} list
#' @param params parameters for the prior distribution. Should match the number and constrains imposed by the distribution name.
#' @param bounds a lower and upper bound for the parameter. This restricts the range in which sampling occurs. If left as NA, no bounds are imposed
#' @return a new \code{prior_dist} object in the arguments specify a valid distribution. Throws an error otherwise.
#' 
#' @export
prior_dist <- function(name, params, bounds = NA){
  if(!is.character(name)){
    stop("prior name must be a string!")
  }
  if(!is.numeric(params)){
    stop("parameters should be numeric!")
  }
  return(validate_prior_dist(new_prior_dist(name, params, bounds)))
}

# Functions below are from https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/qr_regression/stan_utility.R

# Check transitions that ended with a divergence
check_div <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  return(n);
}

# Check transitions that ended prematurely due to maximum tree depth limit
check_exceeded_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  return(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                n, N, max_depth, 100 * n / N))
  if (n > 0)
    return('  Run again with max_depth set to a larger value to avoid saturation')
}

# Checks the potential scale reduction factors
check_rhat <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if ((rhat > 1.1 || is.infinite(rhat)) && !is.nan(rhat)) {
      print(sprintf('Rhat for parameter %s is %s!',
                    rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    return('Rhat looks reasonable for all parameters')
  else
    return('Rhat above 1.1 indicates that the chains very likely have not mixed')
}
