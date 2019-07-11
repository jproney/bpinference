# run this on the DFCI Kraken cluster
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all()
nsims = strtoi(args[1])

e_mat =  matrix(c(2,0),ncol=1)
p_vec = c(1, 1)
z0 = c(1000) # initial population vector
tf = 5 #final simulation timepoint
times = seq(1,tf)

func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5] - c[6]/(1 + exp(c[7]*(x[1] - c[8])))') #logistic function
priors = rep(list(list(name="uniform",params=c(0,2), bounds=c(0,2))),8)
priors[[3]] = list(name="uniform",params=c(5,15), bounds=c(5,15))
priors[[7]] = list(name="uniform",params=c(5,15), bounds=c(5,15))

mod = bp_model(e_mat, p_mat, func_deps, 8, 1)

simulation_params = c(0.15, .2, 12, 0.6, 0.2, .05, 8, .4)

c_mat= matrix(c(0.0,0.2,0.4,0.6,0.8,1.0), ncol=1)

simulation_data = bpsims(mod, simulation_params, z0, times, rep(20,6), c_mat)

dat = stan_data_from_simulation(simulation_data, mod)

generate(mod, priors, "/michorlab/jroney/stanfiles/double_logistic_benchmark.stan")

options(mc.cores = parallel::detectCores())

ranges = matrix(rep(c(0,1), length(simulation_params)),ncol=2,byrow = T)
ranges[3,] = c(5,15)
ranges[7,] = c(5,15)
init = uniform_initialize(ranges, 4)

options(mc.cores = parallel::detectCores())

ranges = matrix(rep(c(0,1), length(simulation_params)),ncol=2,byrow = T)
ranges[3,] = c(5,15)
init = uniform_initialize(ranges, 4)

if(file.exists("/michorlab/jroney/compiles/double_logistic_benchmark.RDS")){
  stan_mod <- readRDS("/michorlab/jroney/compiles/double_logistic_benchmark.RDS")
} else{
  stan_mod <- stan_model(file = "/michorlab/jroney/stanfiles/double_logistic_benchmark.stan")
  saveRDS(stan_mod, "/michorlab/jroney/compiles/double_logistic_benchmark.RDS")
}

fit_data = sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init)
s = extract(fit_data)
ndiv = check_div(fit_data)
poster = data.frame(nsims = nsims, 
                    Theta1_mean = mean(s$Theta1), Theta1_med = median(s$Theta1),  Theta1_95lci = quantile(s$Theta1, .025)[[1]], Theta1_95uci = quantile(s$Theta1, .975)[[1]], Theta1_50lci = quantile(s$Theta1, .25)[[1]], Theta1_50uci = quantile(s$Theta1, .75)[[1]],
                    Theta2_mean = mean(s$Theta2), Theta2_med = median(s$Theta2),  Theta2_95lci = quantile(s$Theta2, .025)[[1]], Theta2_95uci = quantile(s$Theta2, .975)[[1]], Theta2_50lci = quantile(s$Theta2, .25)[[1]], Theta2_50uci = quantile(s$Theta2, .75)[[1]],
                    Theta3_mean = mean(s$Theta3), Theta3_med = median(s$Theta3),  Theta3_95lci = quantile(s$Theta3, .025)[[1]], Theta3_95uci = quantile(s$Theta3, .975)[[1]], Theta3_50lci = quantile(s$Theta3, .25)[[1]], Theta3_50uci = quantile(s$Theta3, .75)[[1]],
                    Theta4_mean = mean(s$Theta4), Theta4_med = median(s$Theta4),  Theta4_95lci = quantile(s$Theta4, .025)[[1]], Theta4_95uci = quantile(s$Theta4, .975)[[1]], Theta4_50lci = quantile(s$Theta4, .25)[[1]], Theta4_50uci = quantile(s$Theta4, .75)[[1]],
                    Theta5_mean = mean(s$Theta5), Theta5_med = median(s$Theta5),  Theta5_95lci = quantile(s$Theta5, .025)[[1]], Theta5_95uci = quantile(s$Theta5, .975)[[1]], Theta5_50lci = quantile(s$Theta5, .25)[[1]], Theta5_50uci = quantile(s$Theta5, .75)[[1]],
                    Theta6_mean = mean(s$Theta6), Theta6_med = median(s$Theta6),  Theta6_95lci = quantile(s$Theta6, .025)[[1]], Theta6_95uci = quantile(s$Theta6, .975)[[1]], Theta6_50lci = quantile(s$Theta6, .25)[[1]], Theta6_50uci = quantile(s$Theta6, .75)[[1]],
                    Theta7_mean = mean(s$Theta7), Theta7_med = median(s$Theta7),  Theta7_95lci = quantile(s$Theta7, .025)[[1]], Theta7_95uci = quantile(s$Theta7, .975)[[1]], Theta7_50lci = quantile(s$Theta7, .25)[[1]], Theta7_50uci = quantile(s$Theta7, .75)[[1]],
                    Theta8_mean = mean(s$Theta8), Theta8_med = median(s$Theta8),  Theta8_95lci = quantile(s$Theta8, .025)[[1]], Theta8_95uci = quantile(s$Theta8, .975)[[1]], Theta8_50lci = quantile(s$Theta8, .25)[[1]], Theta8_50uci = quantile(s$Theta8, .75)[[1]],
                    Divergences = ndiv)


write.table(poster, "/michorlab/jroney/saves/double_logistic_benchmark.csv", sep = ",", col.names = !file.exists("/michorlab/jroney/saves/double_logistic_benchmark.csv"), row.names = FALSE, append = T)
