e_mat <-  matrix(c(2,0),ncol=1)
p_vec <- c(1, 1)
z0 <- c(1000) # initial population vector
tf <- 5 #final simulation timepoint
times <- seq(0,tf)

func_deps <- c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
priors <- rep(list(list(name="normal",params=c(0,.25), bounds=c(0,5))),5)
priors[[3]] <- list(name="uniform",params=c(5,15), bounds=c(5,15))

mod <- bp_model(e_mat, p_vec, func_deps, 5, 1)

simulation_params <- c(0.15, .2, 10, 0.5, 0.10)

c_mat <- matrix(c(0.0,0.2,0.4,0.6,0.8,1.0), ncol=1)

simulation_data <- bpsims(mod, simulation_params, z0, times, rep(20,6), c_mat)

dat <- stan_data_from_simulation(simulation_data, mod, simple_bd = T)

stan_code = generate(mod, priors, simple_bd = T)

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1), length(simulation_params)),ncol=2,byrow = T)
ranges[3,] <- c(5,15)
init <- uniform_initialize(ranges, 4)

stan_mod <- rstan::stan_model(model_code = stan_code)
fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
