e_mat =  matrix(c(2,0),ncol=1)
p_vec = c(1, 1)
z0 = c(1000) # initial population vector
tf = 5 #final simulation timepoint
times = c(1)

func_deps = c('c[1]','c[2]')
priors = rep(list(list(name="normal",params=c(0, .25), bounds=c(0,5))),2)

mod = bp_model(e_mat, p_vec, func_deps, 2, 0)

simulation_params = c(0.25, 0.10)

simulation_dat = bpsims(mod, simulation_params, z0, times, 50)

dat = stan_data_from_simulation(simulation_dat, mod)

generate(mod, priors, "test2.stan")

options(mc.cores = parallel::detectCores())

ranges = matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init = uniform_initialize(ranges, 4)

stan_mod <- stan_model(file = "test2.stan")
fit_data = sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
