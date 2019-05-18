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
    matrix[d,d] b2[d];
    
    //unpack the parameter matrices from theta vector
    E = to_matrix(head(theta, m*d), m, d);
    R = to_matrix(segment(theta,m*d+1, m*d),m,d);
    
    for(i in 1:d){
      b[i] = (R'[i]*E)/sum(col(R,i)); //normalize
    }
    
    
    for(a in 1:d){//the ancestor type we are multiplying by
      for(i in 1:d){
        E_tmp[i] = E[i]*E[i,a];
      }
      
      for(i in 1:d){
        b2[a][i] = (R'[i]*E_tmp)/sum(col(R,i)); //computing the TRANSPOSE of previous convention
      }
    }
    
    for(i in 1:d){
      for(j in 1:d){
        c[i][j] = b2[j][i];
        c[i][j,j] = b2[j][i,j] - b[i,j];//subtract first moment off the diagonal
      }
    }
    
    //unpack the moment matrices from the state vector
    mt = to_matrix(head(state, d*d), d, d, 0);
    dt = to_matrix(segment(state,d*d+1, d*d*d),d,d*d, 1);//read this in row-major order
    
    for(i in 1:d){
      real lamb = sum(col(R,i));
      Beta[i] = to_row_vector(lamb*(mt')*c[i]*mt);
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
  int<lower=0> l; //number of distinct timepoints
  matrix[n,d] pop_vec; //endpoint populations
  matrix[n,d] init_pop; //inital population vectors
  matrix<lower=0, upper=1>[m,d] E; //birth events matrix
  int timesIdx[n]; //index of the duration of each run
  real times[l]; //the actual times
}
transformed data{
  real rdata[0];
  int idata[1];
  real init_state[d*d + d*d*d]; //state is unrolled vector of first + second moments
  real m_init[d,d] = rep_array(0.0,d,d); //first moments
  real d_init[d,d,d] = rep_array(0.0,d,d,d); //second moments
  
  //evaluations of the moments at t=0. Mostly zeros, since type a -> b is impossible in 0 time
  for(i in 1:d){
    m_init[i,i] = 1;
    d_init[i,i,i] = 1;
  }
  
  init_state = append_array(to_array_1d(m_init), to_array_1d(d_init));
  idata[1] = d;
  
}
parameters{
  matrix<lower=0, upper=1>[m,d] R; //birth rate matrix
}
transformed parameters{
  real theta[m*2*d] =  append_array(to_array_1d(E),to_array_1d(R));
}
model{
   real moments[l,d*d + d*d*d]; //raw single-ancestor moments vector evolving over time
   matrix[d,d] m_t[n]; //fisrt moment matrices
   matrix[d,d*d] d_t[n]; //second moments indexing goes (j,k,i) 
   vector[d] Mu_t[n]; //mean vectors for each datapoint
   matrix[d,d] Sigma_t[n]; //covariance matrices for each datapoint
   matrix[d,d] temp; //for copying stuff
  
  //put priors on everything
  for(i in 1:m){
    for(j in 1:d){
      R[i,j] ~ normal(0,.5);
    }
  }
  
  moments = integrate_ode_rk45(moment_ode, init_state, 0, times, theta, rdata, idata);
  
  for(k in 1:n){
    int t = timesIdx[k];
    m_t[n] = to_matrix(head(moments[t],d*d), d,d);
    d_t[n] = to_matrix(segment(moments[t],d*d+1, d*d*d),d,d*d);
    
    //plug in the inital conditions
    Mu_t[n] = (m_t[n]')*to_vector(init_pop[n]);
    temp = to_matrix(init_pop[n]*d_t[n],d,d);
    for(i in 1:d){
      for(j in 1:d){
        Sigma_t[n][i,j] = temp[i,j] - init_pop[n]*(col(m_t[n],i).*col(m_t[n],j)); //subtract to get covariance from 2nd moment
      }
    }

    //and finally...
    pop_vec[n] ~ multi_normal(Mu_t[n], Sigma_t[n]);
  }
}
