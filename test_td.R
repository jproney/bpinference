e_mat <-  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
p_vec <- c(1, 1, 2, 2, 1, 2)
z0 <- c(1000,500) # initial population vector
tf <- 20 #final simulation timepoint
times <- 20
reps <- 10000

simulation_params <- c(0.30, 0.05, 0.03, 0.20, 0.15, 0.10)

r_func <- function(t){sapply(simulation_params, function(x){x*exp(-.2*t)})}

r_func_ub <- function(t){sapply(simulation_params, function(x){x*exp(-.2*t)})}

sim_dat <- replicate(reps, bp_td(e_mat, r_func, r_func_ub, p_vec, z0, times))

calculate_moments_td(e_mat, p_vec, r_func, z0,20)
