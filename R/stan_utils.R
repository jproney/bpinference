#' @import rstan
#' @import dplyr
# wrapper script to aid in branching process inference

# model = the branching process model being estimated pop_vec =
# population vectors at end of each run. Dimensions ndatapts x ntypes init_pop =
# population at beginning of each run.  DImensions ndatapts x ntypes times = the
# length of time elapsing between each initial population and each
# final population. dimensions ndatapts x 1 c_mat = matrix of dependent variables
# that vary from run to run. Dimensions ndatapts x ndep functions = vector of
# function expressions specifying how each rate depends on the
# dependent variables in the c_mat matrix Priors = priors for all of the
# parameters. Should be in the form of a list of vectors Priors should
# be a nparams-dimensional list, and each list entry should have a prior for
# that parameter if c_mat or functions is left as NA, inference is
# performed directly on the rates as parameters

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

# create Stan initialization list by selecting initial values uniformly
# at ranom.  ranges is a nparams x 2 matrix, where nparams is the number of
# parameters. Each row has contains a lower and upper bound for
# initialization.

#' @export
uniform_initialize <- function(ranges, nchains)
{
  return(replicate(nchains, list(list(Theta = apply(ranges, 1, function(s)
  {
    runif(1, s[1], s[2])
  })))))
}

# Helper function for piping simulation output into inference engine

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

# functional_deps is an array of strings encoding functions of variable
# 'x1, x2, ...' and constant parameters 'c1','c2',... priors is list of
# with prior objects for all parameters in order of the rate they
# contribute to example prior object: p =
# list(name='gamma',params=ndep_levels(1,1), init = ndep_levels(0,3), bounds=ndep_levels(0,3))
# filename is name of generated Stan file

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

# takes simple mathematical R expressions and turns them into a
# specific type of Stan code to fill in the template
exp_to_stan <- function(exprn, max_params, max_dep)
{
  if (is.atomic(exprn) || is.name(exprn))
  {
    first <- exprn
  } else
  {
    first <- exprn[[1]]
  }
  
  if (is.atomic(first))
  {
    return(first)
  }
  if (deparse(first) %in% c("+", "-", "/", "*", "^", "exp", "log"))
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
  if (deparse(first) == "(")
  {
    return(exp_to_stan(exprn[[2]], max_params, max_dep))
  }
  if (deparse(first) == "[")
  {
    if (!is.na(stringr::str_extract(deparse(exprn), "^c\\[[1-9]+\\]$")))
    {
      # function parameters
      s <- stringr::str_extract(deparse(exprn), "^c\\[[1-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_params)
      {
        stop(paste("Parameter", s, "goes beyond the number of paramters specified.", 
                   sep = " "))
      }
      return(sprintf("Theta%d", num))
    }
    if (!is.na(stringr::str_extract(deparse(exprn), "^x\\[[1-9]+\\]$")))
    {
      # variables
      s <- stringr::str_extract(deparse(exprn), "^x\\[[1-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_dep)
      {
        stop(paste("Variable", s, "goes beyond the number of dependent variables specified.", 
                   sep = " "))
      }
      return(sprintf("function_var[i,%d]", num))
    }
  } else
  {
    stop("Invalid expression!")
  }
}

check_valid <- function(exprn, max_params, max_dep)
{
  if (is.atomic(exprn) || is.name(exprn))
  {
    first <- exprn
  } else
  {
    first <- exprn[[1]]
  }
  
  if (is.atomic(first))
  {
    return()
  }
  if (deparse(first) %in% c("+", "-", "/", "*", "^", "exp", "log"))
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
  if (deparse(first) == "(")
  {
    check_valid(exprn[[2]], max_params, max_dep)
    return()
  }
  if (deparse(first) == "[")
  {
    if (!is.na(stringr::str_extract(deparse(exprn), "^c\\[[1-9]+\\]$")))
    {
      # function parameters
      s <- stringr::str_extract(deparse(exprn), "^c\\[[1-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_params)
      {
        stop(paste("Parameter", s, "goes beyond the number of paramters specified.", 
                   sep = " "))
      }
      return()
    }
    if (!is.na(stringr::str_extract(deparse(exprn), "^x\\[[1-9]+\\]$")))
    {
      # variables
      s <- stringr::str_extract(deparse(exprn), "^x\\[[1-9]+\\]$")
      num <- strtoi(substr(s, 3, nchar(s) - 1))
      if (num > max_dep)
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

# Functions below are from https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/qr_regression/stan_utility.R

# Check transitions that ended with a divergence
check_div <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  return(n);
}

# Check transitions that ended prematurely due to maximum tree depth limit
check_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                n, N, max_depth, 100 * n / N))
  if (n > 0)
    print('  Run again with max_depth set to a larger value to avoid saturation')
}

# Checks the energy Bayesian fraction of missing information (E-BFMI)
check_energy <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('E-BFMI indicated no pathological behavior')
  else
    print('  E-BFMI below 0.2 indicates you may need to reparameterize your model')
}

# Checks the effective sample size per iteration
check_n_eff <- function(fit) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      print(sprintf('n_eff / iter for parameter %s is %s!',
                    rownames(fit_summary)[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('n_eff / iter looks reasonable for all parameters')
  else
    print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
}

# Checks the potential scale reduction factors
check_rhat <- function(fit) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      print(sprintf('Rhat for parameter %s is %s!',
                    rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('Rhat looks reasonable for all parameters')
  else
    print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
}

check_all_diagnostics <- function(fit) {
  check_n_eff(fit)
  check_rhat(fit)
  check_div(fit)
  check_treedepth(fit)
  check_energy(fit)
}

# Returns parameter arrays separated into divergent and non-divergent transitions
partition_div <- function(fit) {
  nom_params <- extract(fit, permuted=FALSE)
  n_chains <- dim(nom_params)[2]
  params <- as.data.frame(do.call(rbind, lapply(1:n_chains, function(n) nom_params[,n,])))
  
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  params$divergent <- divergent
  
  div_params <- params[params$divergent == 1,]
  nondiv_params <- params[params$divergent == 0,]
  
  return(list(div_params, nondiv_params))
}
