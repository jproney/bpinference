context("test calculation of branching process moments")

test_that("Incorrect inputs to moment calculation function trigger apropriate errors" , {
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  z0 = c(1000) # initial population vector
  simulation_params = c(0.25, 0.10, .5)
  expect_error(calculate_moments(e_mat, p_vec, simulation_params, z0, 5),"Incorrect number of rate parameters provided for model!")

  z0 = c(1000, 500)
  simulation_params = c(0.25, 0.10)
  expect_error(calculate_moments(e_mat, p_vec, simulation_params, z0, 5),"Dimensions of birth events matrix disagree with dimension of initial population vector!")
  
  z0 = c(1000)
  p_vec = c(1,1,1)
  expect_error(calculate_moments(e_mat, p_vec, simulation_params, z0, 5),"Dimensions of birth events matrix disagree with dimension of parent vector!")
  
  p_vec = c(1,2)
  expect_error(calculate_moments(e_mat, p_vec, simulation_params, z0, 5),"Parent vector must have integer entries between 1 and the number of types!")
  
  p_vec = c(0,2)
  expect_error(calculate_moments(e_mat, p_vec, simulation_params, z0, 5),"Parent vector must have integer entries between 1 and the number of types!")
  
  e_mat =  matrix(c(2,0,0,2,0,0,0,0),ncol=2)
  p_vec = c(1, 2, 1.5, 2)
  z0 = c(1000, 500) # initial population vector
  simulation_params = c(0.2, 0.20, .1,.1)  
  expect_error(calculate_moments(e_mat, p_vec, simulation_params, z0, 5),"Parent vector must have integer entries between 1 and the number of types!")
})