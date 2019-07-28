#' model class that can be used for inference and simulation

#' constructor for object of class \code{bp_model}
#' 
#' @param e_mat The matrix of birth events that can occur in the brnaching process. Dimensions \code{nevents} x \code{ntypes}
#' @param p_vec A vector containing the parent type for each of the birth events in \code{e_mat}. Dimensions \code{nevents} x 1
#' @param func_deps A vector of expressions encoding how each rate depends on the parameters of the model and a set of dependent variables.
#' @param nparam The number of parameters in the model
#' @param ndep The number of dependent variables which influence the rates of the model
#' 
#' @return A new branching process model object
new_bp_model <- function(e_mat, p_vec, func_deps, nparam, ndep)
{
  
  return(structure(list(e_mat = e_mat, p_vec = p_vec, func_deps = func_deps, nparams = nparam, 
                        ndep = ndep), class = "bp_model"))
}

#' function to determine whether an object is a valid instance of class \code{bpmodel}
#' 
#' @param bpm  The branching process model object
#' 
#' @return The model passed into the function if it is vaild. Throws an exception otherwise
validate_bp_model <- function(bpm)
{
  if(!is.numeric(bpm$p_vec) || !is.numeric(bpm$e_mat) || !is.numeric(bpm$nparams) ||  !is.numeric(bpm$ndep)){
    stop("e_mat, p_mat, ndep, and nparams must all be numeric!")
  }
  if (!is.vector(bpm$p_vec) || !all(round(bpm$p_vec) == bpm$p_vec) || !all(bpm$p_vec > 0) || !all(bpm$p_vec <= bpm$ntypes))
  {
    stop("p_vec must be a vector of parents for each bith event, where each entry is an positive integer corresponding to the parent type")
  }
  if (is.vector(bpm$e_mat))
  {
    warning("e_mat is a vector, not a matrix. Converting to a one-column matrix")
    bpm$e_mat <- matrix(bpm$e_mat, ncol = 1)
  }
  if (!is.matrix(bpm$e_mat) || !all(round(bpm$e_mat) == bpm$e_mat || !all(bpm$e_mat >=0)))
  {
    stop("e_mat must be a matrix of birth events, where each entry is an integer number of offspring of a given type")
  }
  if (!(is.vector(bpm$func_deps)))
  {
    stop("func_deps must be a vector of strings relating brith rates to dependent variables")
  }
  new_funcs <- sapply(1:length(bpm$func_deps), function(i)
  {
    parse(text = bpm$func_deps[i])[[1]]
  })
  bpm$func_deps <- new_funcs
  if (nrow(bpm$e_mat) != length(bpm$p_vec) || nrow(bpm$e_mat) != length(bpm$func_deps))
  {
    stop("Dimensions of e_mat, p_vec, and func_deps do not agree")
  }
  if (round(bpm$nparams) != bpm$nparams || bpm$nparams <= 0)
  {
    stop("nparams should be an integer number of parameters")
  }
  if (round(bpm$ndep) != bpm$ndep || bpm$ndep < 0)
  {
    stop("ndep should be an integer number of depedent variables")
  }
  for (i in 1:length(bpm$func_deps))
  {
    check_valid(bpm$func_deps[[i]], bpm$nparam, bpm$ndep)
  }
  return(bpm)
}

#' User-friendly function to construct a new valid \code{bpm} object
#' 
#' @param e_mat The matrix of birth events that can occur in the brnaching process. Dimensions \code{nevents} x \code{ntypes}
#' @param p_vec A vector containing the parent type for each of the birth events in \code{e_mat}. Dimensions \code{nevents} x 1
#' @param func_deps A vector of expressions encoding how each rate depends on the parameters of the model and a set of dependent variables.
#' @param nparam The number of parameters in the model
#' @param ndep The number of dependent variables which influence the rates of the model
#' 
#' @return A new branching process model object if the parameters passed are valid. Otherwise an exception is thrown.
#' @export
bp_model <- function(e_mat, p_vec, func_deps, nparam, ndep)
{
  return(validate_bp_model(new_bp_model(e_mat, p_vec, func_deps, nparam, ndep)))
}


#' Convenience constructor for building a simple 1-type birth-death process
#' 
#' @param func_deps A vector of expressions encoding how each rate depends on the parameters of the model and a set of dependent variables.
#' @param nparam The number of parameters in the model
#' @param ndep The number of dependent variables which influence the rates of the model
#' @return A new branching process model object if the parameters passed are valid. Otherwise an exception is thrown.
#' @export
bp_model_simple_birth_death <- function(func_deps, nparam, ndep){
  e_mat = matrix(c(2,0),ncol=1)
  p_vec = c(1,1)
  return(validate_bp_model(new_bp_model(e_mat, p_vec, func_deps, nparam, ndep)))
}
