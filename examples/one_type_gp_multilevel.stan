// Stan file for inferring simple one-type birth death processes wtih functional dependencies. Dispenses with ODE solutions to speed up computation.
data{
  int<lower=1> ndatapts; //number of datapoints
  int<lower=1> ndep_levels; //number of distinct dependent variable levels
  int<lower=1> ndep; //number of dependent variables
  int<lower=1> nparams; //number of parameters
  int<lower=1> nlines; //number of cell lines
  int<lower=1> ndatasets; //number of datasets
  vector[ndatapts] pop_vec; //endpoint populations
  vector[ndatapts] init_pop; //inital population vectors
  int var_idx[ndatapts]; //index of the dependent variable value of each run
  vector[ndep] function_var[ndep_levels]; //the actual dependent variable values
  real times[ndatapts]; //the actual times
  int datasets[ndatapts]; //which datasets point belongs to
  int lines[ndatapts]; //which cell line the point belongs to
}
transformed data{
    vector[ndep_levels] mu = rep_vector(0, ndep_levels);
}
parameters{
    //BEGIN AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
    
	vector[ndep_levels] birth_tilde_global;
	vector[ndep_levels] death_tilde_global;
	vector[ndep_levels] birth_tilde_line[nlines]; 
	vector[ndep_levels] death_tilde_line[nlines]; 
	vector[ndatasets] offsets_b;
	vector[ndatasets] offsets_d;
	
	real<lower=0> rho_b_global;
	real<lower=0> alpha_b_global;

	real<lower=0> rho_d_global;
	real<lower=0> alpha_d_global;

	real<lower=0> rho_b_local;
	real<lower=0> alpha_b_local;

	real<lower=0> rho_d_local;
	real<lower=0> alpha_d_local;

    //END AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
}
transformed parameters{
  matrix[ndep_levels, ndep_levels] K_b_global = cov_exp_quad(function_var, alpha_b_global, rho_b_global) + diag_matrix(rep_vector(1e-8, ndep_levels));
  matrix[ndep_levels, ndep_levels] L_K_b_global = cholesky_decompose(K_b_global);
  matrix[ndep_levels, ndep_levels] K_d_global = cov_exp_quad(function_var, alpha_d_global, rho_d_global) + diag_matrix(rep_vector(1e-8, ndep_levels));
  matrix[ndep_levels, ndep_levels] L_K_d_global = cholesky_decompose(K_d_global);
  vector<lower=0>[ndep_levels] birth_global = exp(L_K_b_global * birth_tilde_global); 
  vector<lower=0>[ndep_levels] death_global = exp(L_K_d_global * death_tilde_global);
  matrix[ndep_levels, ndep_levels] K_b_line[nlines];
  matrix[ndep_levels, ndep_levels] L_K_b_line[nlines];
  matrix[ndep_levels, ndep_levels] K_d_line[nlines];
  matrix[ndep_levels, ndep_levels] L_K_d_line[nlines];
  vector[ndep_levels] birth_line[nlines]; 
  vector[ndep_levels] death_line[nlines]; 
  for(i in 1:nlines){
      K_b_line[i] = cov_exp_quad(function_var, alpha_b_local, rho_b_local) + diag_matrix(rep_vector(1e-8, ndep_levels));
      K_d_line[i] = cov_exp_quad(function_var, alpha_d_local, rho_d_local) + diag_matrix(rep_vector(1e-8, ndep_levels));
      L_K_b_line[i] = cholesky_decompose(K_b_line[i]);
      L_K_d_line[i] = cholesky_decompose(K_d_line[i]);
      birth_line[i] = exp(L_K_b_line[i] * birth_tilde_line[i]); 
      death_line[i] = exp(L_K_d_line[i] * death_tilde_line[i]); 
  }
}
model{
  

  rho_b_global ~ inv_gamma(5, 5);
  alpha_b_global ~ std_normal();

  rho_d_global ~ inv_gamma(5, 5);
  alpha_d_global ~ std_normal();
  
  rho_b_local ~ inv_gamma(5, 5);
  alpha_b_local ~ std_normal();

  rho_d_local ~ inv_gamma(5, 5);
  alpha_d_local ~ std_normal();

  offsets_b ~ normal(0, .1);
  offsets_d ~ normal(0, .1);

  birth_tilde_global ~ std_normal();

  death_tilde_global ~ std_normal();

  for(i in 1:nlines){
    birth_tilde_line[i] ~ std_normal();
    death_tilde_line[i] ~ std_normal(); 
  }
  
  for(k in 1:ndatapts){
    real t = times[k];
    int v = var_idx[k];
    int l = lines[k];
    int ds = datasets[k];
    real b = birth_global[v] + birth_line[l,v] + offsets_b[ds];
    real d = death_global[v] + death_line[l,v] + offsets_d[ds];
    real mu_norm = exp(t*(b-d))*init_pop[k];
    real sigma_norm = sqrt((exp(t*(b-d))*(b*(2*exp(t*(b-d))-1)-d)/(b-d) - exp(2*t*(b-d)))*init_pop[k]); //analytical solution to variance ODE
    
    pop_vec[k] ~ normal(mu_norm, sigma_norm);
  }
}

