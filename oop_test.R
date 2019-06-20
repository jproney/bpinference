d=2
E =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
P = c(1, 1, 2, 2, 1, 2)
Z0 = c(1000,500) # initial population vector
Tf = 5 #final simulation timepoint
times = seq(1,Tf)

func_deps = c('c[1]','c[2]','-c[3]', 'c[4] + x[1]', 'c[5] + exp(log(c[3]))', 'c[6]')
priors = rep(list(list(name="normal",params=c(0,.25), bounds=c(0,5))),6)

mod = bp_model(E, P, func_deps, 6, 1)

generate(mod, priors, "test.stan")

simulation_params = c(0.20, 0.10, 0.05, 0.20, 0.25, 0.20)

bpsims(mod, simulation_params, Z0, times, 50)