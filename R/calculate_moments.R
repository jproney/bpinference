
#' Computes the first and second moments of the the joint population-size distribution of a Markov branching proccess
#'
#' @param e_mat The matrix of birth events that can occur in the brnaching process. Dimensions \code{nevents} x \code{ntypes}
#' @param p_vec A vector containing the parent type for each of the birth events in \code{e_mat}. Dimensions \code{nevents} x 1
#' @param r_vec A vector containing the rate at which each of the birth events in \code{e_mat} occurs. Dimensions \code{nevents} x 1
#' @param z0_vec The initial population vector at time 0. Dimensions \code{ntypes} x 1
#' @return A list of two elements: \code{mu_mat} contains the mean number of each type at time \code{tf}, \code{sigma_mat} is the covariance matrix.
#' 
#' @export
calculate_moments <- function(e_mat,p_vec,r_vec,z0_vec,tf){
  if(!is.numeric(e_mat) || !is.numeric(r_vec) || !is.numeric(z0_vec) || !is.numeric(tf)){
    stop("all parameters should be numeric!")
  }
  if(is.vector(e_mat)){
    warning("e_mat is a vector, not a matrix. Converting to a one-column matrix.")
    e_mat <- matrix(e_mat, ncol = 1)
  }
  if(!is.vector(p_vec) || !is.vector(r_vec) || !is.vector(z0_vec)){
    stop("p_vec, r_vec, and z0_vec should all be vectors!")
  }
  if(tf < 0){
    stop("tf cannot be negative!")
  }
  ntype <- length(z0_vec)
  if(nrow(e_mat) != length(r_vec)){
    stop("Incorrect number of rate parameters provided for model!")
  }
  if(ncol(e_mat) != ntype){
    stop("Dimensions of birth events matrix disagree with dimension of initial population vector!")
  }
  if(nrow(e_mat) != length(p_vec)){
    stop("Dimensions of birth events matrix disagree with dimension of parent vector!")
  }
  if(any(p_vec > ntype | p_vec < 1 | round(p_vec) != p_vec)){
    stop("Parent vector must have integer entries between 1 and the number of types!")
  }
  r_prime_mat <- matrix(rep(0,nrow(e_mat)*ntype), c(nrow(e_mat),ntype))
  r_prime_mat[cbind(1:nrow(e_mat), p_vec)] = r_vec
  
  lamb <- colSums(r_prime_mat)
  b_mat <- t(t(e_mat)%*%r_prime_mat)/lamb
  
  #Diff EQ: m_mat(t) = exp(At)
  a_mat <-  lamb*(b_mat - diag(ntype))
  m_mat <- expm::expm(a_mat*tf)

  c_mat <- array(rep(0,ntype**3), c(ntype, ntype, ntype)); #matrix of second derivatives of offspring PGF
  
  for(i in 1:ntype){
    i_mat <- t(t(e_mat*e_mat[,i])%*%r_prime_mat)/lamb
    i_mat[,i] <- i_mat[,i] - b_mat[,i]
    c_mat[,i,] <- t(i_mat)
    c_mat[i,,] <- t(i_mat)
  }
  
  second_moment_de <- function(t, state, params){
    mt_mat <- expm::expm(a_mat*t)
    beta_mat <- matrix(rep(0,ntype**3), nrow = ntype, ncol = ntype*ntype); #beta_mat array for second moment ODE
    for(i in 1:ntype){
      ai <- lamb[i]
      beta_mat[i,] <- c(ai*t(mt_mat)%*%c_mat[,,i]%*%mt_mat) #vectorized computation of beta_mat
    }
    x <- matrix(state,c(ntype,ntype*ntype))
    x_prime <- c(a_mat%*%x + beta_mat)
    return(list(x_prime))
  }
  
  init_state <- array(rep(0,ntype**3),c(ntype,ntype,ntype)) #inital values of second moments
  for(i in 1:ntype){
    init_state[i,i,i] <- 1;
  }
  init_state <- c(init_state)
  times <- seq(0,tf)
  
  out <- deSolve::ode(y = init_state, times, func = second_moment_de, parms = 0)
  dt_mat <- array(out[nrow(out),-1],c(ntype,ntype*ntype)) #second moment array
  
  mu_mat <- t(m_mat)%*%z0_vec #final mean population matrix
  sigma_mat <- matrix(z0_vec%*%dt_mat,c(ntype,ntype*ntype)) #final population covariance matrix
  
  for(i in 1:ntype){
    for(j in 1:ntype){
      sigma_mat[i,j] <-  sigma_mat[i,j] - z0_vec%*%(m_mat[,i]*m_mat[,j])
    }
  }
  return(list("mu_mat" = mu_mat, "sigma_mat" = sigma_mat, "first_moment_matrix" = m_mat))
}
