#' @import rstan
#' @import dplyr

#' function to generate a stan input list from population data
#' 
#' @param model The branching process model being estimated
#' @param final_pop Population vectors at end of each run. Dimensions ndatapts x ntypes 
#' @param init_pop = population at beginning of each run.  Dimensions ndatapts x ntypes 
#' @param times = the length of time elapsing between each initial population and each final population. dimensions ndatapts x 1 
#' @param c_mat = matrix of dependent variables that vary from run to run. Dimensions ndatapts x ndep 
#'
#' @return A data list to pass to the Stan sampling function 
#' @export
create_stan_data <- function(model, final_pop, init_pop, times, c_mat = NA)
{
  nevents <- nrow(model$e_mat)  #number of events
  ntypes <- ncol(model$e_mat)  #number of types
  ndatapts <- nrow(final_pop)  #number of datapoints
  times_unique <- unique(times)  #distinct timepoints
  ntimes_unique <- length(times_unique)  #number of distinct durations
  times_idx <- match(times, times_unique)  #index of time duration for each datapoint
  nparams <- model$nparams  #total number of parameters
  
  if (model$ndep > 0 && (is.na(c_mat) || ncol(c_mat) < mod$ndep))
  {
    stop("c_mat does not contain enough dependent variables for the model")
  }
  if (model$ndep == 0 && is.na(c_mat))
  {
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
      which.min(abs(rowSums(sweep(c_mat_unique, 2, r))))
    })  #indices of unique dependent variable combinations
  }
  
  
  stan_dat <- list(ntypes = ntypes, nevents = nevents, ndatapts = ndatapts, ntimes_unique = ntimes_unique, ndep_levels = ndep_levels, ndep = ndep, nparams = nparams, e_mat = model$e_mat, 
                   p_vec = model$p_vec, pop_vec = final_pop, init_pop = init_pop, times = array(times_unique, 
                                                                                  1), times_idx = times_idx, function_var = c_mat_unique, var_idx = var_idx)
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
  return(replicate(nchains, list(list(Theta = apply(ranges, 1, function(s)
  {
    runif(1, s[1], s[2])
  })))))
}

#' Helper function for piping simulation output into inference engine
#' @param sim_data simulation data as returned from \code{bpsims}
#' @param model the \code{bpmodel} object the produced the data
#' 
#' @return A data list to pass to the Stan sampling function 
#' @export
stan_data_from_simulation <- function(sim_data, model)
{
  ntypes <- ncol(model$e_mat)
  offset <- model$ndep + (model$ndep > 0)
  for (i in 1:ntypes)
  {
    cellname <- sprintf("t%d_cells", i)
    names(sim_data)[2 + offset + i] <- cellname
    prevname <- paste(cellname, "prev", sep = "_")
    sim_data <- sim_data %>% dplyr::mutate(prev = lag(sim_data[, cellname]))
    names(sim_data)[2 + offset + ntypes + i] <- prevname
  }
  sim_data <- sim_data %>% dplyr::mutate(dtimes = times - lag(times))
  sim_data <- sim_data %>% filter(times != 0)
  init_pop <- as.matrix(sim_data[, (3 + offset + ntypes):(2 + offset + 2 * ntypes)])
  final_pop <- as.matrix(sim_data[, (3 + offset):(2 + offset + ntypes)])
  times <- sim_data$dtimes
  if (model$ndep > 0)
  {
    return(create_stan_data(model, final_pop, init_pop, times, matrix(sim_data[, 
                                                                        4:(3 + model$ndep)], ncol = model$ndep)))
  }
  return(create_stan_data(model, final_pop, init_pop, times))
}

#' Generate a Stan file for inferring the model parameters
#' @param model The \code{bpmodel} object from which to generate the Stan file
#' @param priors A list of prior distributions for each parameters in the mode. 
#' Each entry in the list should be a names list of the form list(name="normal",params=c(0,1),bounds=(-100,100))
#' @param filename The name of the file in which to store the model. Should end in ".stan"
#' @export
generate <- function(model, priors, filename)
{
  
  nfunc <- length(model$func_deps)
  funcs <- rep(0, nfunc)
  
  nprior <- length(priors)
  declStrs <- rep(0, nprior)
  priorStrs <- rep(0, nprior)
  
  for (i in 1:nfunc)
  {
    # convert the expression
    exprn <- model$func_deps[[i]]
    stanstr <- exp_to_stan(exprn, model$nparams, mod$ndep)
    funcs[i] <- sprintf("\t\tr_mat[%d, p_vec[%d]] = %s;\n", i, i, stanstr)
  }
  
  for (p in 1:nprior)
  {
    ind <- which(!is.na(stringr::str_extract(ALLOWED_DISTRIBUTIONS, paste("^", priors[[p]]$name, 
                                                 " ", sep = ""))))
    if (length(ind) > 0)
    {
      nparm <- strtoi(stringr::str_extract(ALLOWED_DISTRIBUTIONS[ind], "[1-9]+"))
      if (nparm != length(priors[[p]]$params))
      {
        stop("Incorrect number of parameters for prior distribution")
      }
      constraints <- stringr::str_extract(ALLOWED_DISTRIBUTIONS[ind], "[_+]+")
      for (i in 1:nparm)
      {
        if (priors[[p]]$params[i] <= 0 && substr(constraints, i, 
                                                 i) == "+")
        {
          stop("Negative parameter not allowed here!")
        }
      }
    } else
    {
      stop("Invalid prior name!")
    }
    priorStrs[p] <- sprintf("\tTheta%d ~ %s(%s);\n", p, priors[[p]]$name, 
                            paste(priors[[p]]$params, collapse = ","))
    if (!any(is.na(priors[[p]]$bounds)))
    {
      declStrs[p] <- sprintf("\treal<lower=%d, upper=%d> Theta%d;\n", 
                             priors[[p]]$bounds[1], priors[[p]]$bounds[2], p)
    }
  }
  write(sprintf(STAN_TEMPLATE, paste(declStrs, collapse = ""), paste(funcs, 
                                                                collapse = ""), paste(priorStrs, collapse = "")), filename)
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
