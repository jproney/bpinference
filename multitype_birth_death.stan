functions{
  //one giant ODE system for solving all the moements at once
  //theta = unrolled birth matrix + death matrix
  //state - unrolled vector of first moment matrix + second moment matrix
  real[] moment_ode(real t, real[] state, real[] theta, real[] rdata, int[] idata) {// rdata param not used
    int d = idata[1];
    int m = idata[2];
    matrix[m,d] E;             //birth events
    matrix[m,d] R;             //birth rates
    matrix[d,d] A;             //1st moment ode coefficients
    matrix[d,d] b;             //individual offspring mean matrix
    matrix[d,d] c[d];          //array of matrices of second derivates of offspring PGFs
    matrix[d,d*d] Beta;        //2nd moment ode coefficients
    matrix[d,d] mt;            //mean matrix at time t
    matrix[d,d] mt_prime;      //first moment derivatives
    matrix[d,d*d] dt;          //second moments at time t. Each col is an (j,k) covariance pair. Row i is ancestor
    matrix[d,d*d] dt_prime;    //second moment derivatives. Each col is an (j,k) covariance pair. Row i is ancestor
    matrix[m,d] E_tmp;
    matrix[d,d] b2;
    real lambda[d];
    
    //unpack the parameter matrices from theta vector
    E = to_matrix(head(theta, m*d), m, d);
    R = to_matrix(segment(theta,m*d+1, m*d),m,d);
    
    for(i in 1:d){
      lambda[i] = sum(col(R,i));
      b[i] = (R'[i]*E)/lambda[i]; //normalize
      A[i] = b[i]*lambda[i];
      A[i,i] = (b[i,i] - 1)*lambda[i];
    }
    
    for(a in 1:d){//the ancestor type we are multiplying by
      for(i in 1:m){
        E_tmp[i] = E[i]*E[i,a];
      }
      
      for(i in 1:d){
        b2[i] = (R'[i]*E_tmp)/lambda[i]; //computing the TRANSPOSE of previous convention
        c[i][a] = b2[i];
        c[i][a,a] = b2[i,a] - b[i,a];
      }
    }
    
    
    //unpack the moment matrices from the state vector
    mt = to_matrix(head(state, d*d), d, d, 0);
    dt = to_matrix(segment(state,d*d+1, d*d*d),d,d*d, 1);//read this in row-major order
    
    for(i in 1:d){
      Beta[i] = to_row_vector(lambda[i]*(mt')*c[i]*mt);
    }

    //first moment ODE
    mt_prime = A*mt;

    //second moment ODE
    dt_prime = A*dt + Beta;
    
    //collapse derivatives into 1d array and return
    return append_array(to_array_1d(mt_prime),to_array_1d(dt_prime));
  }
}
data{
  int<lower=0> d; //number of types
  int<lower=0> m; //number of events
  int<lower=0> n; //number of datapoints
  int<lower=0> l; //number of distinct timepoints.
  int<lower=0> c; //number of distinct dependent variable levels
  matrix[n,d] pop_vec; //endpoint populations
  matrix[n,d] init_pop; //inital population vectors
  matrix<lower=0>[m,d] E; //birth events matrix
  int<lower=0> P[m]; //parents for each birth event
  vector<lower=0>[m] Pri_mu; //prior means for each rate
  vector<lower=0>[m] Pri_sig; //prior variances for each rate
  int<lower=0> func_type[m]; //functional dependence. 1 = constant (1 param), 2 = linear (2 params), 3 = quadratic (3 params), 4 = logistic (4 params)
  int var_idx[n]; //index of the dependent variable value of each run
  real function_var[c]; //the actual dependent variable values
  int times_idx[n]; //index of the duration of each run
  real times[l]; //the actual times
}
transformed data{
  int total_parms = sum(func_type);
  int event_idx[m]; //vector of indices for which parameters pretain to which birth event
  int count = 1;
  real rdata[0];
  int idata[2];
  real init_state[d*d + d*d*d]; //state is unrolled vector of first + second moments
  real m_init[d,d] = rep_array(0.0,d,d); //first moments
  real d_init[d,d,d] = rep_array(0.0,d,d,d); //second moments
  
  for(i in 1:m){
    event_idx[i] = count;
    count = count + func_type[i];
  }
  
  //evaluations of the moments at t=0. Mostly zeros, since type a -> b is impossible in 0 time
  for(i in 1:d){
    m_init[i,i] = 1;
    d_init[i,i,i] = 1;
  }
  
  init_state = append_array(to_array_1d(m_init), to_array_1d(d_init));
  idata[1] = d;
  idata[2] = m;
  
}
parameters{
  vector<lower=0, upper=2>[total_parms] R; //parameters of birth/death rate functions
}
transformed parameters{
  matrix[m,d] R_prime;
  real theta[c,m*2*d];
  for(i in 1:c){
    R_prime = rep_matrix(0, m,d);
    for(k in 1:m){
      if(func_type[k] == 1){
          R_prime[k, P[k]] = R[event_idx[k]]; //constant parameter
      }
      else if(func_type[k] == 2){
          R_prime[k,P[k]] = R[event_idx[k]] + R[event_idx[k]+1]*function_var[l]; //linear response variable
      }
      else if(func_type[k] == 3){
          R_prime[k,P[k]] = R[event_idx[k]] + R[event_idx[k]+1]*function_var[l] + R[event_idx[k]+1]*(function_var[l]^2); //quadratic response variable
      }
      else if(func_type[k] == 4){
          R_prime[k,P[k]] = R[event_idx[k]] + (R[event_idx[k]+1] - R[event_idx[k]])/(1 + (function_var[l] / R[event_idx[k]+2])^R[event_idx[k]+3]); //logistic response variable
      }
    }
    theta[c] =  append_array(to_array_1d(E),to_array_1d(R_prime));
  }
}
model{
   real moments[c,l,d*d + d*d*d]; //raw single-ancestor moments vector evolving over time
   matrix[d,d] m_t; //fisrt moment matrices
   matrix[d,d*d] d_t; //second moments indexing goes (j,k,i) 
   vector[d] Mu_t; //mean vectors for each datapoint
   matrix[d,d] Sigma_t; //covariance matrices for each datapoint
   matrix[d,d] temp; //for copying stuff
  
  
  //put priors on everything
  R ~ normal(Pri_mu, Pri_sig);
  
  for(i in 1:c){
    moments[c] = integrate_ode_rk45(moment_ode, init_state, 0, times, theta[i], rdata, idata);
  }
  
  for(k in 1:n){
    int t = times_idx[k];
    int v = var_idx[k];
    m_t = to_matrix(head(moments[v,t],d*d), d,d);
    d_t = to_matrix(segment(moments[v,t],d*d+1, d*d*d),d,d*d);
    
    //plug in the inital conditions
    Mu_t = (m_t')*to_vector(init_pop[k]);
    temp = to_matrix(init_pop[k]*d_t,d,d);
    for(i in 1:d){
      for(j in 1:d){
        Sigma_t[i,j] = temp[i,j] - init_pop[k]*(col(m_t,i).*col(m_t,j)); //subtract to get covariance from 2nd moment
      }
    }

    //and finally...
    pop_vec[k] ~ multi_normal(Mu_t, Sigma_t);
  }
}
