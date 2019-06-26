# run this on the DFCI Kraken cluster
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all()
nsims = strtoi(args[1])
e_mat =  matrix(c(2,0),ncol=1)
p_vec = c(1, 1)
z0 = c(1000) # initial population vector
tf = 5 #final simulation timepoint
times = seq(1,tf)

func_deps = c('c[1]','c[2]')
priors = rep(list(list(name="uniform",params=c(0, 2), bounds=c(0,2))),2)

mod = bp_model(e_mat, p_vec, func_deps, 2, 0)

simulation_params = c(0.25, 0.10)

simulation_dat = bpsims(mod, simulation_params, z0, times, nsims)

dat = stan_data_from_simulation(simulation_dat, mod)

generate(mod, priors, "/michorlab/jroney/stanfiles/one_type_benchmark_unif.stan")

options(mc.cores = parallel::detectCores())

ranges = matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init = uniform_initialize(ranges, 4)

if(file.exists("/michorlab/jroney/compiles/one_type_benchmark_unif.RDS")){
  stan_mod <- readRDS("/michorlab/jroney/compiles/one_type_benchmark_unif.RDS")
} else{
  stan_mod <- stan_model(file = "/michorlab/jroney/stanfiles/one_type_benchmark_unif.stan")
  saveRDS(stan_mod, "/michorlab/jroney/compiles/one_type_benchmark_unif.RDS")
}

fit_data = sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
s = extract(fit_data)
poster = data.frame(nsims = nsims, Theta1_mean = mean(s$Theta1),Theta2_mean = mean(s$Theta2), Theta1_med = median(s$Theta1), Theta2_med = median(s$Theta2),
                    Theta1_95lci = quantile(s$Theta1, .025)[[1]], Theta2_95lci = quantile(s$Theta2, .025)[[1]], 
                    Theta1_95uci = quantile(s$Theta1, .975)[[1]], Theta2_95uci = quantile(s$Theta2, .975)[[1]])

write.table(poster, "/michorlab/jroney/saves/one_type_benchmark_unif.csv", sep = ",", col.names = !file.exists("/michorlab/jroney/saves/one_type_benchmark_unif.csv"), row.names = FALSE, append = T)
