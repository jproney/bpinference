context("test simulation of multitype branching processes")

test_that("Simulation runs without error and returns data of proper size -- no dependencies case",{
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  z0 = c(1000) # initial population vector
  tf = 5 #final simulation timepoint
  times = seq(0,tf)
  reps=5
  
  func_deps = c('c[1]','c[2]')
  mod = bp_model(e_mat, p_vec, func_deps, 2, 0)
  simulation_params = c(0.25, 0.10)
  simulation_dat = bpsims(mod, simulation_params, z0, times, reps)
  expect_equal(simulation_dat$pop[1], z0)
  expect_equal(nrow(simulation_dat), length(times)*reps)
  
  e_mat =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
  p_vec = c(1, 1, 2, 2, 1, 2)
  z0 = c(1000,500) # initial population vector
  tf = 5 #final simulation timepoint
  times = seq(0,tf)
  reps=5
  
  func_deps = c('c[1]','c[2]','c[3]', 'c[4]', 'c[5]', 'c[6]')
  mod = bp_model(e_mat, p_vec, func_deps, 6, 0)
  simulation_params = c(0.20, 0.10, 0.05, 0.20, 0.25, 0.20)
  simulation_dat = bpsims(mod, simulation_params, z0, times, reps)
  expect_equal(simulation_dat[1,3], z0[1])
  expect_equal(simulation_dat[1,4], z0[2])
  expect_equal(nrow(simulation_dat), length(times)*reps)
  
})

test_that("Simulation runs without error and returns data of proper size -- cases with functional dependencies",{
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  z0 = c(1000) # initial population vector
  tf = 5 #final simulation timepoint
  times = seq(0,tf)
  reps = replicate(6,1+rpois(1,4))
  
  func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
  mod = bp_model(e_mat, p_vec, func_deps, 5, 1)
  simulation_params = c(0.15, .2, 10, 0.5, 0.10)
  c_mat = matrix(c(0.0,0.2,0.4,0.6,0.8,1.0), ncol=1)
  simulation_dat = bpsims(mod, simulation_params, z0, times, reps, c_mat)
  expect_equal(nrow(simulation_dat), length(times)*sum(reps))
  expect_equal(simulation_dat$pop[1], z0)
  expect_equal(simulation_dat$dep[1], c_mat[1])
  expect_equal(simulation_dat$dep[nrow(simulation_dat)], c_mat[length(c_mat)])
  
  e_mat =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
  p_vec = c(1, 1, 2, 2, 1, 2)
  z0 = c(1000,1000) # initial population vector
  tf = 5 #final simulation timepoint
  times = seq(0,tf)
  reps = replicate(6,1+rpois(1,4))
  
  func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]','c[6]', 'c[7] + c[8]/(1 + exp(c[9]*(x[2] - c[10])))', 'c[11]', 'c[12]')
  mod = bp_model(e_mat, p_vec, func_deps, 12, 2)
  simulation_params = c(0.15, .2, 10, 0.5, 0.05, 0.05, 0.15, .2, 10, 0.5, 0.10, 0.10)
  c_mat = matrix(rep(c(0.0,0.2,0.4,0.6,0.8,1.0),2),ncol=2,nrow=6)
  simulation_dat = bpsims(mod, simulation_params, z0, times, reps, c_mat)
  expect_equal(nrow(simulation_dat), length(times)*sum(reps))
  expect_equal(simulation_dat[1,4], c_mat[1,1])
  expect_equal(simulation_dat[1,5], c_mat[1,2])
  expect_equal(simulation_dat[1,6], z0[1])
  expect_equal(simulation_dat[1,7], z0[2])
  expect_equal(simulation_dat[nrow(simulation_dat),4], c_mat[nrow(c_mat),1])
  expect_equal(simulation_dat[nrow(simulation_dat),5], c_mat[nrow(c_mat),2])
})