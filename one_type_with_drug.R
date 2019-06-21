E =  matrix(c(2,0),ncol=1)
P = c(1, 1)
Z0 = c(1000) # initial population vector
Tf = 5 #final simulation timepoint
times = seq(1,Tf)

func_deps = c('c[1]','c[2]')
priors = rep(list(list(name="normal",params=c(0,.25), bounds=c(0,5))),2)

mod = bp_model(E, P, func_deps, 2, 0)

simulation_params = c(0.35, 0.10)

X = bpsims(mod, simulation_params, Z0, times, 50)

dat = stan_data_from_simulation(X, mod)

generate(mod, priors, "test2.stan")

options(mc.cores = parallel::detectCores())

ranges = matrix(rep(c(0,1),nrow(E)),nrow(E),2,byrow = T)
init = uniform_initialize(ranges, 4)

library(rstan)
stan_mod <- stan_model(file = "test2.stan")
fit.data = sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
