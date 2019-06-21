
#Computes the first and second moment of the population vector at time Tf for a 
#branching process with parameters E, P, R and initial population 

#' @export
calculate_moments <- function(E,P,R,Z0,Tf){
  
  d <- length(Z0)
  R_prime <- matrix(rep(0,nrow(E)*d), c(nrow(E),d))
  R_prime[cbind(1:nrow(E), P)] = R
  
  lamb <- colSums(R_prime)
  b <- t(t(E)%*%R_prime)/lamb
  
  #Diff EQ: M(t) = exp(At)
  A <-  lamb*(b - diag(d))
  M <- expm::expm(A*Tf)
  
  C <- array(rep(0,d**3), c(d, d, d)); #matrix of second derivatives of offspring PGF
  
  for(i in 1:d){
    i_mat <- t(t(E*E[,i])%*%R_prime)/lamb
    i_mat[,i] <- i_mat[,i] - b[,i]
    C[,i,] <- t(i_mat)
    C[i,,] <- t(i_mat)
  }
  
  second_moment_de <- function(t, state, params){
    mt <- expm::expm(A*t)
    Beta <- matrix(rep(0,d**3), nrow = d, ncol = d*d); #Beta array for second moment ODE
    for(i in 1:d){
      ai <- lamb[i]
      Beta[i,] <- c(ai*t(mt)%*%C[,,i]%*%mt) #vectorized computation of Beta
    }
    x <- matrix(state,c(d,d*d))
    x_prime <- c(A%*%x + Beta)
    return(list(x_prime))
  }
  
  init_state <- array(rep(0,d**3),c(d,d,d)) #inital values of second moments
  for(i in 1:d){
    init_state[i,i,i] <- 1;
  }
  init_state <- c(init_state)
  times <- seq(0,Tf)
  
  out <- deSolve::ode(y = init_state, times, func = second_moment_de, parms = 0)
  dt <- array(out[nrow(out),-1],c(d,d*d)) #second moment array
  
  Mu <- t(M)%*%Z0 #final mean population matrix
  Sigma <- matrix(Z0%*%dt,c(d,d*d)) #final population covariance matrix
  
  for(i in 1:d){
    for(j in 1:d){
      Sigma[i,j] <-  Sigma[i,j] - Z0%*%(M[,i]*M[,j])
    }
  }
  return(list("mu" = Mu, "sigma" = Sigma))
}
