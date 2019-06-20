d=2
E =  matrix(c(2,0),2,1)
P = c(1,1)
Z0 = c(1000) # initial population vector
Tf = 5 #final simulation timepoint

f1 = function(x){.15 + .2*exp(-2*x)}
f2 = function(x){.2}
fs = c(f1,f2)

C = matrix(c(0.0,0.2,0.4,0.6,0.8,1.0),6,1)

times = seq(1,Tf)
X = bpsims_modeling(E, P, Z0, C, fs, times = times, reps = rep(50,6))
X <- X %>% group_by(rep)
names(X)[4] <- "cells"
names(X)[5] <- "drug_level"
X <- X %>% group_by(rep) %>% mutate(cells_prev = lag(cells), dtimes = times - lag(times))

# Remove NA values
X <- X %>% filter(!is.na(cells_prev))

pop_vec = X$cells
init_pop = X$cells_prev

f = c(quote(c1 + c2*exp(-c3*x1)), quote(c4))
prior = list(list(name = "uniform", params=c(0,1)), list(name = "uniform", params=c(0,1)), list(name = "lognormal", params=c(0,1)), list(name = "uniform", params=c(0,1)))
ranges = matrix(c(0,1,0,1,0,3,0,1),4,2,byrow = T)

out = create_stan_data(E = E, P= P, final_pop = matrix(pop_vec,ncol=1), init_pop = matrix(init_pop, ncol=1), times = X$dtimes, priors = prior, C = matrix(X$drug_level, ncol=1), functions = f)
init = uniform_initialize(ranges, 4)

options(mc.cores = parallel::detectCores())
fit.data = sampling(out$model, data = out$data, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
