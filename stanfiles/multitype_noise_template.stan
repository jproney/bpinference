functions{
  //one giant ODE system for solving all the moements at once
  //theta = unrolled birth matrix + death matrix
  //state - unrolled vector of first moment matrix + second moment matrix
  real[] moment_ode(real t, real[] state, real[] theta, real[] rdata, int[] idata) {// rdata param not used
    int ntypes = idata[1];
    int nevents = idata[2];
    matrix[nevents,ntypes] e_mat;             //birth events
    matrix[nevents,ntypes] r_mat;             //birth rates
    matrix[ntypes,ntypes] a_mat;             //1st moment ode coefficients
    matrix[ntypes,ntypes] b_mat;             //individual offspring mean matrix
    matrix[ntypes,ntypes] c_mat[ntypes];          //array of matrices of second derivates of offspring PGFs
    matrix[ntypes,ntypes*ntypes] beta_mat;        //2nd moment ode coefficients
    matrix[ntypes,ntypes] mt;            //mean matrix at time t
    matrix[ntypes,ntypes] mt_prime;      //first moment derivatives
    matrix[ntypes,ntypes*ntypes] dt;          //second moments at time t. Each col is an (j,k) covariance pair. Row i is ancestor
    matrix[ntypes,ntypes*ntypes] dt_prime;    //second moment derivatives. Each col is an (j,k) covariance pair. Row i is ancestor
    matrix[nevents,ntypes] e_mat_tmp;
    matrix[ntypes,ntypes] b2;
    real lambda[ntypes];
    
    //unpack the parameter matrices from theta vector
    e_mat = to_matrix(head(theta, nevents*ntypes), nevents, ntypes);
    r_mat = to_matrix(segment(theta,nevents*ntypes+1, nevents*ntypes),nevents,ntypes);

    for(i in 1:ntypes){
      lambda[i] = sum(col(r_mat,i));
      b_mat[i] = (r_mat'[i]*e_mat);
      a_mat[i] = b_mat[i];
      a_mat[i,i] = b_mat[i,i] - lambda[i];
    }
    
    for(a in 1:ntypes){//the ancestor type we are multiplying by
      for(i in 1:nevents){
        e_mat_tmp[i] = e_mat[i]*e_mat[i,a];
      }
      
      for(i in 1:ntypes){
        b2[i] = (r_mat'[i]*e_mat_tmp); //computing the TRANSPOSE of previous convention
        c_mat[i][a] = b2[i];
        c_mat[i][a,a] = b2[i,a] - b_mat[i,a];
      }
    }
    
    //unpack the moment matrices from the state vector
    mt = to_matrix(head(state, ntypes*ntypes), ntypes, ntypes);
    dt = to_matrix(segment(state,ntypes*ntypes+1, ntypes*ntypes*ntypes),ntypes,ntypes*ntypes);

    for(i in 1:ntypes){
      beta_mat[i] = to_row_vector((mt')*c_mat[i]*mt);
    }

    //first moment ODE
    mt_prime = a_mat*mt;
    
    //second moment ODE
    dt_prime = a_mat*dt + beta_mat;

    //collapse derivatives into 1d array and return
    return append_array(to_array_1d(mt_prime),to_array_1d(dt_prime));
  }
}
data{
  int<lower=0> ntypes; //number of types
  int<lower=0> nevents; //number of events
  int<lower=0> ndatapts; //number of datapoints
  int<lower=0> ntimes_unique; //number of distinct timepoints.
  int<lower=0> ndep_levels; //number of distinct dependent variable levels
  int<lower=0> ndep; //number of dependent variables
  int<lower=0> nparams; //number of parameters
  matrix[ndatapts,ntypes] pop_vec; //endpoint populations
  matrix[ndatapts,ntypes] init_pop; //inital population vectors
  matrix<lower=0>[nevents,ntypes] e_mat; //birth events matrix
  int<lower=0> p_vec[nevents]; //parents for each birth event
  int var_idx[ndatapts]; //index of the dependent variable value of each run
  real function_var[ndep_levels,ndep]; //the actual dependent variable values
  int times_idx[ndatapts]; //index of the duration of each run
  real times[ntimes_unique]; //the actual times
}
transformed data{
  real rdata[0];
  int idata[2];
  real init_state[ntypes*ntypes + ntypes*ntypes*ntypes]; //state is unrolled vector of first + second moments
  real m_init[ntypes,ntypes] = rep_array(0.0,ntypes,ntypes); //first moments
  real d_init[ntypes,ntypes,ntypes] = rep_array(0.0,ntypes,ntypes,ntypes); //second moments
  
  //evaluations of the moments at t=0. Mostly zeros, since type a -> b_mat is impossible in 0 time
  for(i in 1:ntypes){
    m_init[i,i] = 1;
    d_init[i,i,i] = 1;
  }
  
  init_state = append_array(to_array_1d(m_init), to_array_1d(d_init));
  idata[1] = ntypes;
  idata[2] = nevents;
  
}
parameters{
    //BEGIN AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
    
%s
    //END AUTO-GENERATED CODE FOR PARAMETER DECLARATIONS
    real<lower=0> sigma_obs;
}
transformed parameters{
  matrix[nevents,ntypes] r_mat;
  real theta[ndep_levels,nevents*2*ntypes];
  for(i in 1:ndep_levels){
    r_mat = rep_matrix(0, nevents,ntypes);
    
    //BEGIN AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
    
%s
    //END AUTO-GENERATED CODE FOR FUNCTIONAL DEPENDENCIES
    theta[i] =  append_array(to_array_1d(e_mat),to_array_1d(r_mat));
  }
}
model{
  real moments[ndep_levels,ntimes_unique,ntypes*ntypes + ntypes*ntypes*ntypes]; //raw single-ancestor moments vector evolving over time
  matrix[ntypes,ntypes] m_t; //fisrt moment matrices
  matrix[ntypes,ntypes*ntypes] d_t; //second moments indexing goes (j,k,i) 
  vector[ntypes] mu_t; //mean vectors for each datapoint
  matrix[ntypes,ntypes] sigma_t; //covariance matrices for each datapoint
  matrix[ntypes,ntypes] temp; //for copying stuff
  
  
  //BEGIN AUTO-GENERATED CODE FOR PRIORS

%s
  //END AUTO-GENERATED CODE FOR PRIORS
  sigma_obs ~ normal(0,.1);
  
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
    pop_vec[k] ~ multi_normal(mu_t, sigma_t + sigma_obs*diag_matrix(mu_t));
  }
}

