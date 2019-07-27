// Stan file for inferring simple one-type birth death processes wtih functional dependencies. Dispenses with ODE solutions to speed up computation.
data{
  int<lower=0> ndatapts; //number of datapoints
  int<lower=0> ndep_levels; //number of distinct dependent variable levels
  int<lower=0> ndep; //number of dependent variables
  int<lower=0> nparams; //number of parameters
  vector[ndatapts] pop_vec; //endpoint populations
  vector[ndatapts] init_pop; //inital population vectors
  int var_idx[ndatapts]; //index of the dependent variable value of each run
  real function_var[ndep_levels,ndep]; //the actual dependent variable values
  real times[ndatapts]; //the actual times
}
parameters{
    //BEGIN AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
    
%s
    //END AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
}
transformed parameters{
  matrix[2,ndep_levels] r_mat;
  for(i in 1:ndep_levels){
    //BEGIN AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
    
%s
    //END AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
  }
}
model{
  
  
  //BEGIN AUTO-GENERATED CODE FOR PRIORS
  
%s
  //END AUTO-GENERATED CODE FOR PRIORS
  
  
  for(k in 1:ndatapts){
    real t = times[k];
    int v = var_idx[k];
    real b = r_mat[1,v];
    real d = r_mat[2,v];
    real mu = exp(t*(b-d))*init_pop[k];
    real sigma = sqrt((exp(t*(b-d))*(b*(2*exp(t*(b-d))-1)-d)/(b-d) - exp(2*t*(b-d)))*init_pop[k]); //analytical solution to variance ODE
    
    pop_vec[k] ~ normal(mu, sigma);
  }
}
