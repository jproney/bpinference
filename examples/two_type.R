e_mat <-  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
p_vec <- c(1, 1, 2, 2, 1, 2)
z0 <- c(1000,500) # initial population vector
tf <- 5 #final simulation timepoint
times <- seq(0,tf)

func_deps <- c('c[1]','c[2]','c[3]', 'c[4]', 'c[5]', 'c[6]')
priors <- rep(list(list(name="normal",params=c(0,.25), bounds=c(0,5))),6)

mod <- bp_model(e_mat, p_vec, func_deps, 6, 0)

simulation_params <- c(0.20, 0.10, 0.05, 0.20, 0.25, 0.20)

simulation_dat <- bpsims(mod, simulation_params, z0, times, 50)

dat <- stan_data_from_simulation(simulation_dat, mod)

stan_code = generate(mod, priors, loo = T)

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init <- uniform_initialize(ranges, 1)

stan_mod <- rstan::stan_model(model_code = stan_code)
fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.8), chains = 1, refresh = 1, init =init, iter=3000,warmup=1000)
