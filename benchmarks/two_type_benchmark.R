# run this on the DFCI Kraken cluster
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all()
nsims = strtoi(args[1])

e_mat =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
p_vec = c(1, 1, 2, 2, 1, 2)
z0 = c(1000,500) # initial population vector
tf = 5 #final simulation timepoint
times = seq(1,tf)

func_deps = c('c[1]','c[2]','c[3]', 'c[4]', 'c[5]', 'c[6]')
priors = rep(list(list(name="uniform",params=c(0,2), bounds=c(0,2))),6)

mod = bp_model(e_mat, p_vec, func_deps, 6, 0)

simulation_params = c(0.20, 0.10, 0.05, 0.20, 0.25, 0.20)

simulation_dat = bpsims(mod, simulation_params, z0, times, nsims)

dat = stan_data_from_simulation(simulation_dat, mod)

generate(mod, priors, "/michorlab/jroney/stanfiles/two_type_benchmark.stan")

ranges = matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init = uniform_initialize(ranges, 4)


if(file.exists("/michorlab/jroney/compiles/two_type_benchmark.RDS")){
  stan_mod <- readRDS("/michorlab/jroney/compiles/two_type_benchmark.RDS")
} else{
  stan_mod <- stan_model(file = "/michorlab/jroney/stanfiles/two_type_benchmark.stan")
  saveRDS(stan_mod, "/michorlab/jroney/compiles/two_type_benchmark.RDS")
}

options(mc.cores = parallel::detectCores())

fit_data = sampling(stan_mod, data = dat, control = list(adapt_delta = 0.8), chains = 4, refresh = 1, init =init)
s = extract(fit_data)
poster = data.frame(nsims = nsims, 
                    Theta1_mean = mean(s$Theta1), Theta1_med = median(s$Theta1),  Theta1_95lci = quantile(s$Theta1, .025)[[1]], Theta1_95uci = quantile(s$Theta1, .975)[[1]], 
                    Theta2_mean = mean(s$Theta2), Theta2_med = median(s$Theta2),  Theta2_95lci = quantile(s$Theta2, .025)[[1]], Theta2_95uci = quantile(s$Theta2, .975)[[1]], 
                    Theta3_mean = mean(s$Theta3), Theta3_med = median(s$Theta3),  Theta3_95lci = quantile(s$Theta3, .025)[[1]], Theta3_95uci = quantile(s$Theta3, .975)[[1]], 
                    Theta4_mean = mean(s$Theta4), Theta4_med = median(s$Theta4),  Theta4_95lci = quantile(s$Theta4, .025)[[1]], Theta4_95uci = quantile(s$Theta4, .975)[[1]], 
                    Theta5_mean = mean(s$Theta5), Theta5_med = median(s$Theta5),  Theta5_95lci = quantile(s$Theta5, .025)[[1]], Theta5_95uci = quantile(s$Theta5, .975)[[1]], 
                    Theta6_mean = mean(s$Theta6), Theta6_med = median(s$Theta6),  Theta6_95lci = quantile(s$Theta6, .025)[[1]], Theta6_95uci = quantile(s$Theta6, .975)[[1]])
                    

write.table(poster, "/michorlab/jroney/saves/two_type_benchmark.csv", sep = ",", col.names = !file.exists("/michorlab/jroney/saves/two_type_benchmark.csv"), row.names = FALSE, append = T)
