d=2
E =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
P = c(1, 1, 2, 2, 1, 2)
R = c(.2, .1, .05, .2, .25, .2)
Z0 = c(1000,500) # initial population vector
Tf = 5 #final simulation timepoint

times = seq(1,Tf)
X = bpsims(E, R, P, Z0,times = times, reps = 50)
X <- X %>% group_by(rep)
names(X)[3:4] <- c("t1_cells", "t2_cells")
X <- X %>% group_by(rep) %>% mutate(t1_cells_prev = lag(t1_cells), t2_cells_prev = lag(t2_cells),
                                    t1_cells_0 = first(t1_cells), t2_cells_0 = first(t2_cells),
                                    dtimes = times - lag(times))

# Remove NA values
X <- X %>% filter(!is.na(t1_cells_prev))

pop_vec = cbind(X$t1_cells,X$t2_cells)
init_pop = cbind(X$t1_cells_prev,X$t2_cells_prev)
model_dtimes = sort(unique(times[-1] - times[-length(times)]))

mom = calculate_moments(E,P,R,Z0,Tf)

Pri_mu = rep(0, nrow(E))
Pri_sig = rep(.25, nrow(E))
scale = 1/max(Pri_sig)

library(rstan)
options(mc.cores = parallel::detectCores())
stan_dat <- list(d = d, m = nrow(E), n = nrow(pop_vec), l=1,  pop_vec = pop_vec, init_pop = init_pop, E = E, P = P, 
                 times = as.array(model_dtimes)/scale, times_idx = as.numeric(factor(X$dtimes)), Pri_mu = Pri_mu*scale, Pri_sig = Pri_sig*scale, 
                      c = 1, func_type = rep(1,nrow(E)), var_idx = rep(1,nrow(pop_vec)), function_var = array(0,1))

fit <- stan_model(file = "multitype_birth_death.stan")
fit.data <- sampling(fit, data = stan_dat, control = list(adapt_delta = 0.8), chains = 4, refresh = 1)
s = extract(fit.data)

t = format(Sys.time(), "%a_%b_%d_%H:%M:%S_%Y")
saveRDS(fit.data, paste("saves/",t, "_",toString(d), "type_samples.rds", sep=""))
saveRDS(R, paste("saves/",t, "_", toString(d), "type_truth.rds", sep=""))
saveRDS(stan_dat, paste("saves/",t, "_", toString(d), "type_data.rds", sep=""))