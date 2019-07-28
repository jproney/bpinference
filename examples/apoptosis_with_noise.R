devtools::load_all()
detect_p = .9

e_mat <-  rbind(c(2,0),c(0,1),c(0,0))
p_vec <- c(1, 1, 2)
z0 <- c(1000,200) # initial population vector
times <- seq(0,5)

func_deps <- c('c[1]','c[2]', 'c[3]')
priors <- rep(list(list(name="normal",params=c(0, .25), bounds=c(0,2))),3)

mod <- bp_model(e_mat, p_vec, func_deps, 3, 0)

simulation_params <- c(0.08, 0.05, .1)

simulation_dat <- bpsims(mod, simulation_params, z0, times, 50)

simulation_dat$X1 = sapply(simulation_dat$X1, function(x){rbinom(1,x, detect_p)}) #add measurement noise
simulation_dat$X2 = sapply(simulation_dat$X2, function(x){rbinom(1,x, detect_p)}) #add measurement noise


dat <- stan_data_from_simulation(simulation_dat, mod)


stan_code <- generate(mod, priors, noise_model = T)

stan_mod <- stan_model(model_code = stan_code)
options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init <- uniform_initialize(ranges, 4)

fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, iter = 3000, warmup = 1000, init=init)
