context("test correctness of bpmodel class")

test_that("Incorrect objects are not instantiated",{
  e_mat =  matrix(c(2,0),ncol=1)
  p_mat = c(1, 1)
  func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5] - c[6]/(1 + exp(c[7]*(x[1] - c[11])))') #logistic function
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Parameter c\\[11\\] goes beyond the number of parameters specified.")
  
  e_mat =  matrix(c(2,3,0),ncol=1)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]','c[3]')
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Dimensions of e_mat, p_vec, and func_deps do not agree")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1, 1)
  func_deps = c('c[1]','c[2]')
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Dimensions of e_mat, p_vec, and func_deps do not agree")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, -1)
  func_deps = c('c[1]','c[2]')
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "p_vec must be a vector of parents for each bith event, where each entry is an positive integer corresponding to the parent type")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]', 'c[3]')
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1), "Dimensions of e_mat, p_vec, and func_deps do not agree")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]*x[3]') 
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1),  "Variable x\\[3\\] goes beyond the number of dependent variables specified.")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1.5)
  func_deps = c('c[1]','c[2]*x[3]') 
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 4),  "p_vec must be a vector of parents for each bith event, where each entry is an positive integer corresponding to the parent type")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]*x[3]') 
  expect_error(bp_model(e_mat, p_mat, func_deps, -8, 1),  "nparams should be an integer number of parameters")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('hi','mom') 
  expect_error(bp_model(e_mat, p_mat, func_deps, 8, 1),  "Invalid expression: hi")
  
  e_mat =  matrix(c(2,0,1,1),ncol=2,nrow=2)
  p_mat = c(1, 1)
  func_deps = c('c[1]','c[2]*x[3]') 
  expect_error(bp_model("all", p_mat, func_deps, 2, 1),  "e_mat, p_mat, ndep, and nparams must all be numeric!")
  expect_error(bp_model(e_mat, "your", func_deps, 2, 1),  "e_mat, p_mat, ndep, and nparams must all be numeric!")
  expect_error(bp_model(e_mat, p_mat, func_deps, "base", 1),  "e_mat, p_mat, ndep, and nparams must all be numeric!")
  expect_error(bp_model(e_mat, p_mat, func_deps, 2, "arebelongtous"),  "e_mat, p_mat, ndep, and nparams must all be numeric!")
  expect_error(bp_model(e_mat, list(1,1), func_deps, 2, 1),  "e_mat, p_mat, ndep, and nparams must all be numeric!")
  expect_error(bp_model(e_mat, p_mat, print, 2, 1),  "func_deps must be a vector of strings relating brith rates to dependent variables")
})

test_that("Valid objects are properly instantiated",{
  e_mat =  matrix(c(2,0),ncol=1)
  p_mat = c(1, 1)
  func_deps = c('c[1] + c[2]/(1 + exp(c[3] * (x[1] - c[4])))','c[5] - c[6]/(1 + exp(c[7] * (x[1] - c[8])))') #logistic function
  expect_silent(bp_model(e_mat, p_mat, func_deps, 8, 1))
  mod = bp_model(e_mat, p_mat, func_deps, 8, 1)
  expect_equal(mod$e_mat, e_mat)
  expect_equal(mod$p_vec, p_mat)
  expect_equal(deparse(mod$func_deps[[1]]), func_deps[1])
  expect_equal(deparse(mod$func_deps[[2]]), func_deps[2])
})

