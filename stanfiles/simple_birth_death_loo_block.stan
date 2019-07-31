generated quantities{
  vector[ndatapts] log_lik;
  for(k in 1:ndatapts){
    real t = times[k];
    int v = var_idx[k];
    real b = r_mat[1,v];
    real d = r_mat[2,v];
    real mu = exp(t*(b-d))*init_pop[k];
    real sigma = sqrt((exp(t*(b-d))*(b*(2*exp(t*(b-d))-1)-d)/(b-d) - exp(2*t*(b-d)))*init_pop[k]); //analytical solution to variance ODE
    
    log_lik[k] =  normal_lpdf(pop_vec[k] | mu, sigma);
  }
}