# model class that can be used for inference and simulation

new_bp_model <- function(e_mat, p_vec, func_deps, nparam, ndep)
{
  
  return(structure(list(e_mat = e_mat, p_vec = p_vec, func_deps = func_deps, nparams = nparam, 
                        ndep = ndep), class = "bp_model"))
}

validate_bp_model <- function(bpm)
{
  data <- unclass(bpm)
  if (!is.vector(data$p_vec) || !all(round(data$p_vec) == data$p_vec) || !all(data$p_vec > 
                                                                  0))
  {
    stop("p_vec must be a vector of parents for each bith event, where each entry is an positive integer corresponding to the parent type")
  }
  if (!is.matrix(data$e_mat) || !all(round(data$e_mat) == data$e_mat || !all(data$e_mat >= 
                                                                 0)))
  {
    stop("e_mat must be a matrix of birth events, where each entry is an integer number of offspring of a given type")
  }
  if (!(is.vector(data$func_deps)))
  {
    stop("func_deps must be a vector of language objects relating brith rates to dependent variables")
  }
  
  if (nrow(data$e_mat) != length(data$p_vec) || nrow(data$e_mat) != length(data$func_deps))
  {
    stop("Dimensions of e_mat, p_vec, and func_deps do not agree")
  }
  if (round(data$nparams) != data$nparams || data$nparams <= 0)
  {
    stop("nparams should be an integer number of parameters")
  }
  if (round(data$ndep) != data$ndep || data$ndep < 0)
  {
    stop("ndep should be an integer number of depedent variables")
  }
  for (i in 1:length(func_deps))
  {
    check_valid(data$func_deps[[i]], data$nparam, data$ndep)
  }
  return(bpm)
}

#' @export
bp_model <- function(e_mat, p_vec, func_deps, nparam, ndep)
{
  if (is.vector(e_mat))
  {
    e_mat <- matrix(e_mat, ncol = 1)
  }
  new_funcs <- sapply(1:length(func_deps), function(i)
  {
    parse(text = func_deps[i])[[1]]
  })
  return(validate_bp_model(new_bp_model(e_mat, p_vec, new_funcs, nparam, ndep)))
}

