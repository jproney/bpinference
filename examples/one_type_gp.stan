// Stan file for inferring simple one-type birth death processes wtih functional dependencies. Dispenses with ODE solutions to speed up computation.
data{
  int<lower=1> ndatapts; //number of datapoints
  int<lower=1> ndep_levels; //number of distinct dependent variable levels
  int<lower=1> ndep; //number of dependent variables
  int<lower=1> nparams; //number of parameters
  vector[ndatapts] pop_vec; //endpoint populations
  vector[ndatapts] init_pop; //inital population vectors
  int var_idx[ndatapts]; //index of the dependent variable value of each run
  vector[ndep] function_var[ndep_levels]; //the actual dependent variable values
  real times[ndatapts]; //the actual times
}
transformed data{
    vector[ndep_levels] mu = rep_vector(0, ndep_levels);
}
parameters{
    //BEGIN AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
    
	vector<lower=0, upper=5>[ndep_levels] birth;
	real<lower=0, upper=5> death;
	
	real<lower=0> rho;
	real<lower=0> alpha;
	real<lower=0> sigma;


    //END AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
}
transformed parameters{
  matrix[2,ndep_levels] r_mat;
  for(i in 1:ndep_levels){
    //BEGIN AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
    
		r_mat[1, i] = birth[i];
		r_mat[2, i] = death;


    //END AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
  }
}
model{
  
  matrix[ndep_levels, ndep_levels] L_K;
  matrix[ndep_levels, ndep_levels] K = cov_exp_quad(function_var, alpha, rho);
  real sq_sigma = square(sigma);

  // diagonal elements
  for (n in 1:ndep_levels){
    K[n, n] = K[n, n] + sq_sigma;
  }
  
  //BEGIN AUTO-GENERATED CODE FOR PRIORS

  L_K = cholesky_decompose(K);

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();

  birth ~ multi_normal_cholesky(mu, L_K);

  //END AUTO-GENERATED CODE FOR PRIORS
  
  
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

