#args <- commandArgs(trailingOnly = TRUE)
devtools::load_all()
detect_p = as.numeric(args[1])
detect_p = .8

e_mat <-  rbind(c(2,0),c(0,1),c(0,0))
p_vec <- c(1, 1, 2)
z0 <- c(1000,200) # initial population vector
times <- seq(0,5)

func_deps <- c('c[1]','c[2]', 'c[3]')
priors <- rep(list(list(name="normal",params=c(0, .25), bounds=c(0,2))),3)

mod <- bp_model(e_mat, p_vec, func_deps, 3, 0)

simulation_params <- c(0.08, 0.05, .1)

simulation_dat <- bpsims(mod, simulation_params, z0, times, 25)

simulation_dat$X1 = sapply(simulation_dat$X1, function(x){rbinom(1,x, detect_p)}) #add measurement noise
simulation_dat$X2 = sapply(simulation_dat$X2, function(x){rbinom(1,x, detect_p)}) #add measurement noise


dat <- stan_data_from_simulation(simulation_dat, mod)


#if(file.exists("/michorlab/jroney/compiles/one_type_noise_benchmark.RDS")){
#  stan_mod <- readRDS("/michorlab/jroney/compiles/one_type_noise_benchmark.RDS")
#} else{
#  stan_mod <- stan_model(model_code = stan_code)
#  saveRDS(stan_mod, "/michorlab/jroney/compiles/one_type_noise_benchmark.RDS")
#}
stan_mod <- stan_model(file = "with_apop_noise_model.stan")
options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init <- uniform_initialize(ranges, 4)

fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, iter = 3000, warmup = 1000, init=init)
s = extract(fit_data)

#poster = data.frame(noise = detect_p, Theta1_mean = mean(s$Theta1),Theta2_mean = mean(s$Theta2), Theta1_med = median(s$Theta1), Theta2_med = median(s$Theta2),
#                    Theta1_95lci = quantile(s$Theta1, .025)[[1]], Theta2_95lci = quantile(s$Theta2, .025)[[1]], 
#                    Theta1_95uci = quantile(s$Theta1, .975)[[1]], Theta2_95uci = quantile(s$Theta2, .975)[[1]])

#write.table(poster, "/michorlab/jroney/saves/one_type_noise_benchmark.csv", sep = ",", col.names = !file.exists("/michorlab/jroney/saves/one_type_noise_benchmark.csv"), row.names = FALSE, append = T)

