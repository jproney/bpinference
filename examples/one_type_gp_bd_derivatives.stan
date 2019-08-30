// Stan file for inferring simple one-type birth death processes wtih functional dependencies. Dispenses with ODE solutions to speed up computation.
data{
  int<lower=1> ndatapts; //number of datapoints
  int<lower=1> ntotal; //total number of datapoints, observed and predicted
  int<lower=1> nobs; //number of observed datapoints
  vector[ndatapts] pop_vec; //endpoint populations
  vector[ndatapts] init_pop; //inital population vectors
  int var_idx[ndatapts]; //index of the dependent variable value of each run
  real x_vals[ntotal]; //the actual independent variable values
  real times[ndatapts]; //the actual times
}
transformed data{
    vector[ntotal] mu = rep_vector(0, ntotal);
    int slope_sign[ntotal] = rep_array(0,ntotal);
}
parameters{
    //BEGIN AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
    
	vector[2*ntotal] birth_tilde;
	vector[2*ntotal] death_tilde;
	
	real<lower=0> rho_b;
	real<lower=0> alpha_b;

	real<lower=0> rho_d;
	real<lower=0> alpha_d;


    //END AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
}
transformed parameters{
  matrix[2*ntotal, 2*ntotal] K_b; 
  matrix[2*ntotal, 2*ntotal] L_K_b;
  matrix[2*ntotal, 2*ntotal] K_d;
  matrix[2*ntotal, 2*ntotal] L_K_d;
  vector[ntotal] birth; 
  vector[ntotal] death;
  vector[ntotal] birth_prime; 
  vector[ntotal] death_prime;
  
  matrix[2,nobs] r_mat;
  
  for(i in 1:ntotal){ //Covariance between funtion values and derivative values
    for(j in 1:ntotal){
      K_b[i,j] = alpha_b^2*exp(-(x_vals[i] - x_vals[j])^2/(2*rho_b^2));
    }
  }

  for(i in 1:ntotal){ //Covariance between funtion values and derivative values
    for(j in 1:ntotal){
      K_b[ntotal + i,j] = - alpha_b^2/rho_b^2 * (x_vals[i] - x_vals[j]) * exp(-(x_vals[i] - x_vals[j])^2/(2*rho_b^2));
      K_b[j,ntotal + i] = - alpha_b^2/rho_b^2 * (x_vals[i] - x_vals[j]) * exp(-(x_vals[i] - x_vals[j])^2/(2*rho_b^2));
    }
  }
  
  for(i in 1:ntotal){
    for(j in 1:ntotal){
        K_b[ntotal + i,ntotal + j] = alpha_b^2/rho_b^2 *(1 - 1/rho_b^2*(x_vals[i] - x_vals[j])^2)* exp(-(x_vals[i] - x_vals[j])^2/(2*rho_b^2));
    }
  }
  
  L_K_b = cholesky_decompose(K_b + diag_matrix(rep_vector(1e-8, 2*ntotal)));

  for(i in 1:ntotal){ //Covariance between funtion values and derivative values
    for(j in 1:ntotal){
      K_d[i,j] = alpha_d^2*exp(-(x_vals[i] - x_vals[j])^2/(2*rho_d^2));
    }
  }

  for(i in 1:ntotal){ //Covariance between funtion values and derivative values
    for(j in 1:ntotal){
      K_d[ntotal + i,j] = - alpha_d^2/rho_d^2 * (x_vals[i] - x_vals[j]) * exp(-(x_vals[i] - x_vals[j])^2/(2*rho_d^2));
      K_d[j,ntotal + i] = - alpha_d^2/rho_d^2 * (x_vals[i] - x_vals[j]) * exp(-(x_vals[i] - x_vals[j])^2/(2*rho_d^2));
    }
  }
  
  for(i in 1:ntotal){
    for(j in 1:ntotal){
        K_d[ntotal + i, ntotal + j] = alpha_d^2/rho_d^2 *(1 - 1/rho_d^2*(x_vals[i] - x_vals[j])^2) * exp(-(x_vals[i] - x_vals[j])^2/(2*rho_d^2)); 
    }
  }
  
  L_K_d = cholesky_decompose(K_d+ diag_matrix(rep_vector(1e-8, 2*ntotal)));
  
  birth = head(L_K_b * birth_tilde, ntotal); 
  birth_prime = tail(L_K_b * birth_tilde, ntotal); 
  death = head(L_K_d * death_tilde, ntotal); 
  death_prime = tail(L_K_d * death_tilde, ntotal); 
  
  for(i in 1:nobs){
    //BEGIN AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
    
		r_mat[1, i] = exp(birth[i]);
		r_mat[2, i] = exp(death[i]);


    //END AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
  }
}
model{
  

  rho_b ~ inv_gamma(5, 5);
  alpha_b ~ std_normal();

  birth_tilde ~ std_normal();
  
  rho_d ~ inv_gamma(5, 5);
  alpha_d ~ std_normal();

  death_tilde ~ std_normal();

  slope_sign ~ bernoulli_logit(50*(birth_prime .*exp(birth) - death_prime .*exp(death))); //hyperparameter for punishing positive slopes
  
  for(k in 1:ndatapts){
    real t = times[k];
    int v = var_idx[k];
    real b = r_mat[1,v];
    real d = r_mat[2,v];
    real mu_norm = exp(t*(b-d))*init_pop[k];
    real sigma_norm = sqrt((exp(t*(b-d))*(b*(2*exp(t*(b-d))-1)-d)/(b-d) - exp(2*t*(b-d)))*init_pop[k]); //analytical solution to variance ODE
    
    pop_vec[k] ~ normal(mu_norm, sigma_norm);
  }
}

