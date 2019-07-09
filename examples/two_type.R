e_mat =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
p_vec = c(1, 1, 2, 2, 1, 2)
z0 = c(1000,500) # initial population vector
tf = 5 #final simulation timepoint
times = seq(0,1)

func_deps = c('c[1]','c[2]','c[3]', 'c[4]', 'c[5]', 'c[6]')
priors = rep(list(list(name="normal",params=c(0,.25), bounds=c(0,5))),6)

mod = bp_model(e_mat, p_vec, func_deps, 6, 0)

simulation_params = c(0.25, 0.10, 0.05, 0.20, 0.30, 0.15)

simulation_dat = bpsims(mod, simulation_params, z0, times, 50)

dat = stan_data_from_simulation(simulation_dat, mod)

#reg = solve(t(dat$init_pop)%*%dat$init_pop)%*%t(dat$init_pop)%*%dat$pop_vec

#z_bar = apply(dat$init_pop,2,mean)

#mu_hat = z_bar%*%reg

#sample_var = t(apply(dat$init_pop%*%reg - dat$pop_vec, 1 , function(s){s%*%t(s)}))

#mom2 = solve(t(dat$init_pop)%*%dat$init_pop)%*%t(dat$init_pop)%*%sample_var

#sig_hat = matrix(mom2[1,]*z_bar[1] + mom2[2,]*z_bar[2],2,2,byrow = T)

mom = calculate_moments(e_mat, p_vec, simulation_params, z0, 1)

generate(mod, priors, "test.stan")

options(mc.cores = parallel::detectCores())

ranges = matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init = uniform_initialize(ranges, 4)

stan_mod <- stan_model(file = "test.stan")
fit_data = sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
