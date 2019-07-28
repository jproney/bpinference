e_mat <-  matrix(c(2,0),ncol=1)
p_vec <- c(1, 1)
z0 <- c(1000) # initial population vector
times <- seq(1,20)

func_deps <- c('c[1]','c[2]')
priors <- rep(list(list(name="uniform",params=c(0, 1), bounds=c(0,2))),2)

mod <- bp_model_simple_birth_death(func_deps, 2, 0)

simulation_params <- c(0.1, 0.05)

simulation_dat <- bpsims(mod, simulation_params, z0, times, 25)

dat <- stan_data_from_simulation(simulation_dat, mod, simple_bd = T)

stan_code = generate(mod, priors, simple_bd = T)

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init <- uniform_initialize(ranges, 4)

stan_mod <- rstan::stan_model(model_code = stan_code)
fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, iter = 3000, warmup = 1000)
