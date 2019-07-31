generated quantities {
  vector[ndatapts] log_lik;

  real moments[ndep_levels,ntimes_unique,ntypes*ntypes + ntypes*ntypes*ntypes]; //raw single-ancestor moments vector evolving over time
  matrix[ntypes,ntypes] m_t; //fisrt moment matrices
  matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i) 
  vector[ntypes] mu_t; //mean vectors for each datapoint
  matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
  matrix[ntypes,ntypes] temp; //for copying stuff
  
  for(i in 1:ndep_levels){
    moments[i] = integrate_ode_rk45(moment_ode, init_state, 0, times, theta[i], rdata, idata);
  }
  
  for(k in 1:ndatapts){
    int t = times_idx[k];
    int v = var_idx[k];
    m_t = to_matrix(head(moments[v,t],ntypes*ntypes), ntypes,ntypes);
    d_t = to_matrix(segment(moments[v,t],ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);
    
    //plug in the inital conditions
    mu_t = (m_t')*to_vector(init_pop[k]);
    temp = to_matrix(init_pop[k]*d_t,ntypes,ntypes);
    for(i in 1:ntypes){
      for(j in 1:ntypes){
        sigma_t[i,j] = temp[i,j] - init_pop[k]*(col(m_t,i).*col(m_t,j)); //subtract to get covariance from 2nd moment
      }
    }
    
    //and finally...
    log_lik[k] = multi_normal_lpdf(pop_vec[k] | mu_t, sigma_t);
  }
}

