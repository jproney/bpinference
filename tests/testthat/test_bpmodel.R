context("test correctness of bpmodel class")

test_that("Incorrect objects are not instantiated",{
  e_mat =  matrix(c(2,0),ncol=1)
  p_mat = c(1, 1)
  func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5] - c[6]/(1 + exp(c[7]*(x[1] - c[11])))') #logistic function
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Parameter c\\[11\\] goes beyond the number of parameters specified.")
  
  e_mat =  matrix(c(2,3,0),ncol=1)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]','c[3]') #logistic function
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Dimensions of e_mat, p_vec, and func_deps do not agree")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1, 1)
  func_deps = c('c[1]','c[2]') #logistic function
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Dimensions of e_mat, p_vec, and func_deps do not agree")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]', 'c[3]') #logistic function
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Dimensions of e_mat, p_vec, and func_deps do not agree")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]*x[3]') #logistic function
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1),  "Variable x\\[3\\] goes beyond the number of dependent variables specified.")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]*x[3]') #logistic function
  expect_error(bp_model(e_mat, p_mat, func_deps, -8, 1),  "nparams should be an integer number of parameters")
})
