context("test conversion of real data to stan-friendly format")

test_that("Incorrect inputs trigger apropriate errors" , {
  e_mat <-  matrix(c(2,0),ncol=1)
  p_vec <- c(1, 1)
  z0 <- c(1000) # initial population vector
  times <- seq(1,5)
  
  func_deps <- c('c[1]','c[2]')
  priors <- rep(list(list(name="uniform",params=c(0, 1), bounds=c(0,2))),2)
  
  mod <- bp_model_simple_birth_death(func_deps, 2, 0)
  
  simulation_params <- c(0.1, 0.05)
  simulation_dat <- bpsims(mod, simulation_params, z0, times, 25)
  dat <- stan_data_from_simulation(simulation_dat, mod, simple_bd = T)
  
  expect_error(create_stan_data(mod, final_pop = dat$pop_vec, init_pop = dat$init_pop, times = dat$times[-1], simple_bd = T),"times, initial populations, and final populations must all have same number of rows!")
  expect_error(create_stan_data(mod, final_pop = dat$pop_vec, init_pop = dat$init_pop[-1], times = dat$times, simple_bd = T),"times, initial populations, and final populations must all have same number of rows!")
  expect_error(create_stan_data(mod, final_pop = dat$pop_vec[-1], init_pop = dat$init_pop, times = dat$times, simple_bd = T),"times, initial populations, and final populations must all have same number of rows!")
  expect_warning(create_stan_data(mod, final_pop = dat$pop_vec, init_pop = dat$init_pop, times = dat$times, simple_bd = T),"final_pop is a vector, not a matrix. Converting to a one-column matrix")
  
})
