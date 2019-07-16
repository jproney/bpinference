
#Computes the first and second moment of the population vector at time tf for a_mat 
#branching process with parameters e_mat, p_vec, r_vec and initial population 

#' @export
calculate_moments <- function(e_mat,p_vec,r_vec,z0_vec,tf){
  
  ntype <- length(z0_vec)
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
      print(c_mat[,,i])
      beta_mat[i,] <- c(ai*t(mt_mat)%*%c_mat[,,i]%*%mt_mat) #vectorized computation of beta_mat
    }
    x <- matrix(state,c(ntype,ntype*ntype))
    x_prime <- c(a_mat%*%x + beta_mat)
    print(beta_mat)
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
  return(list("mu_mat" = mu_mat, "sigma_mat" = sigma_mat))
}
