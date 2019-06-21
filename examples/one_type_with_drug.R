E =  matrix(c(2,0),ncol=1)
P = c(1, 1)
Z0 = c(1000) # initial population vector
Tf = 5 #final simulation timepoint
times = seq(1,Tf)

func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
priors = rep(list(list(name="normal",params=c(0,.25), bounds=c(0,5))),5)
priors[[3]] = list(name="uniform",params=c(5,15), bounds=c(5,15))

mod = bp_model(E, P, func_deps, 5, 1)

simulation_params = c(0.15, .2, 10, 0.5, 0.10)

C = matrix(c(0.0,0.2,0.4,0.6,0.8,1.0), ncol=1)

X = bpsims(mod, simulation_params, Z0, times, rep(50,6), C)

dat = stan_data_from_simulation(X, mod)

generate(mod, priors, "test3.stan")

options(mc.cores = parallel::detectCores())

ranges = matrix(rep(c(0,1), length(simulation_params)),ncol=2,byrow = T)
ranges[3,] = c(5,15)
init = uniform_initialize(ranges, 4)

library(rstan)
stan_mod <- stan_model(file = "test3.stan")
fit.data = sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
