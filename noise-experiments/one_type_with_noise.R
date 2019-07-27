args <- commandArgs(trailingOnly = TRUE)
devtools::load_all()
#detect_p = as.numeric(args[1])
detect_p = .8

e_mat <-  matrix(c(2,0),ncol=1)
p_vec <- c(1, 1)
z0 <- c(1000) # initial population vector
times <- seq(1,20)

func_deps <- c('c[1]','c[2]')
priors <- rep(list(list(name="uniform",params=c(0, 1), bounds=c(0,2))),2)

mod <- bp_model_simple_birth_death(func_deps, 2, 0)

simulation_params <- c(0.08, 0.05)

simulation_dat <- bpsims(mod, simulation_params, z0, times, 25)

simulation_dat$pop = sapply(simulation_dat$pop, function(x){rbinom(1,x, detect_p)}) #add measurement noise


dat <- stan_data_from_simulation(simulation_dat, mod, simple_bd = T)

stan_code = generate(mod, priors, "one_type_noise.stan", simple_bd=T)

if(file.exists("/michorlab/jroney/compiles/one_type_noise_benchmark.RDS")){
  stan_mod <- readRDS("/michorlab/jroney/compiles/one_type_noise_benchmark.RDS")
} else{
  stan_mod <- stan_model(model_code = stan_code)
  saveRDS(stan_mod, "/michorlab/jroney/compiles/one_type_noise_benchmark.RDS")
}


stan_mod <- stan_model(file = "one_type_noise.stan")


options(mc.cores = parallel::detectCores())

fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, iter = 3000, warmup = 1000)
s = extract(fit_data)

poster = data.frame(noise = detect_p, Theta1_mean = mean(s$Theta1),Theta2_mean = mean(s$Theta2), Theta1_med = median(s$Theta1), Theta2_med = median(s$Theta2),
                    Theta1_95lci = quantile(s$Theta1, .025)[[1]], Theta2_95lci = quantile(s$Theta2, .025)[[1]], 
                    Theta1_95uci = quantile(s$Theta1, .975)[[1]], Theta2_95uci = quantile(s$Theta2, .975)[[1]])

write.table(poster, "/michorlab/jroney/saves/one_type_noise_benchmark.csv", sep = ",", col.names = !file.exists("/michorlab/jroney/saves/one_type_noise_benchmark.csv"), row.names = FALSE, append = T)

