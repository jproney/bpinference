d=2
E =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
P = c(1, 1, 2, 2, 1, 2)
Z0 = c(1000,500) # initial population vector
Tf = 5 #final simulation timepoint
times = seq(1,Tf)

func_deps = c('c[1]','c[2]','c[3]', 'c[4]', 'c[5] + c[2]', 'c[6] + exp(c[1])')

mod = bp_model(E, P, func_deps, 6, 0)

simulation_params =  c(0.20, 0.10, 0.05, 0.20, 0.25, 0.20)

bpsims_modeling(mod, simulation_params, Z0, times, 50)