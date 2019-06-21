# model class that can be used for inference and simulation

new_bp_model <- function(E, P, func_deps, nParam, nDep)
{
  
  return(structure(list(E = E, P = P, func_deps = func_deps, nParams = nParam, 
                        nDep = nDep), class = "bp_model"))
}

validate_bp_model <- function(bpm)
{
  data <- unclass(bpm)
  if (!is.vector(data$P) || !all(round(data$P) == data$P) || !all(data$P > 
                                                                  0))
  {
    stop("P must be a vector of parents for each bith event, where each entry is an positive integer corresponding to the parent type")
  }
  if (!is.matrix(data$E) || !all(round(data$E) == data$E || !all(data$E >= 
                                                                 0)))
  {
    stop("E must be a matrix of birth events, where each entry is an integer number of offspring of a given type")
  }
  if (!(is.vector(data$func_deps)))
  {
    stop("func_deps must be a vector of language objects relating brith rates to dependent variables")
  }
  
  if (nrow(data$E) != length(data$P) || nrow(data$E) != length(data$func_deps))
  {
    stop("Dimensions of E, P, and func_deps do not agree")
  }
  if (round(data$nParams) != data$nParams || data$nParams <= 0)
  {
    stop("nParams should be an integer number of parameters")
  }
  if (round(data$nDep) != data$nDep || data$nDep < 0)
  {
    stop("nDep should be an integer number of depedent variables")
  }
  for (i in 1:length(func_deps))
  {
    check_valid(data$func_deps[[i]], data$nParam, data$nDep)
  }
  return(bpm)
}

#' @export
bp_model <- function(E, P, func_deps, nParam, nDep)
{
  if (is.vector(E))
  {
    E <- matrix(E, ncol = 1)
  }
  new_funcs <- sapply(1:length(func_deps), function(i)
  {
    parse(text = func_deps[i])[[1]]
  })
  return(validate_bp_model(new_bp_model(E, P, new_funcs, nParam, nDep)))
}

