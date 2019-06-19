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

mom = calculate_moments(E,P,R,Z0,Tf)

prior = rep(list(list(name = "exponential", params=c(2))),nrow(E))
ranges = matrix(rep(c(0,2),nrow(E)),nrow(E),2,byrow = T)

out = create_stan_data(E = E, P= P, final_pop = pop_vec, init_pop = init_pop, times = X$dtimes, priors = prior)
init = uniform_initialize(ranges, 4)

#options(mc.cores = parallel::detectCores())
#fit.data = sampling(out$model, data = out$data, control = list(adapt_delta = 0.8), chains = 4, refresh = 1, init =init)