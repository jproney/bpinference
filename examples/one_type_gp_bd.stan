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
    
	vector[ndep_levels] birth_tilde;
	vector[ndep_levels] death_tilde;
	
	real<lower=0> rho_b;
	real<lower=0> alpha_b;

	real<lower=0> rho_d;
	real<lower=0> alpha_d;


    //END AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
}
transformed parameters{
  matrix[ndep_levels, ndep_levels] K_b = cov_exp_quad(function_var, alpha_b, rho_b) + diag_matrix(rep_vector(1e-8, ndep_levels));
  matrix[ndep_levels, ndep_levels] L_K_b = cholesky_decompose(K_b);
  matrix[ndep_levels, ndep_levels] K_d = cov_exp_quad(function_var, alpha_d, rho_d) + diag_matrix(rep_vector(1e-8, ndep_levels));
  matrix[ndep_levels, ndep_levels] L_K_d = cholesky_decompose(K_d);
  vector<lower=0>[ndep_levels] birth = exp(L_K_b * birth_tilde); 
  vector<lower=0>[ndep_levels] death = exp(L_K_d * death_tilde); 
  
  matrix[2,ndep_levels] r_mat;
  for(i in 1:ndep_levels){
    //BEGIN AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
    
		r_mat[1, i] = birth[i];
		r_mat[2, i] = death[i];


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
