context("test that full stan models are properly generated")

test_that("incorrect models are not generated", {
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  func_deps = c('c[1]','c[2]')
  mod = bp_model(e_mat, p_vec, func_deps, 2, 0)
  priors = rep(list(prior_dist(name="uniform",params=c(0, 1), bounds=c(0,2))),1)
  expect_error(generate(mod,priors, "test_gen.stan"), "incorrect number of priors for model!")
  
  priors = rep(list(prior_dist(name="uniform",params=c(0, 1), bounds=c(0,2))),3)
  expect_error(generate(mod,priors, "test_gen.stan"), "incorrect number of priors for model!")
  
  priors = rep(list(prior_dist(name="uniform",params=c(0, 1), bounds=c(0,2))),2)
  generate(mod,priors, "test_gen.stan")
  out = readChar("test_gen.stan", file.info("test_gen.stan")$size)
  expect_equal(out, REFERENCE_MODEL_1)
  file.remove("test_gen.stan")
})
  
test_that("generated models match known reference models", {
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  func_deps = c('c[1]','c[2]')
  mod = bp_model(e_mat, p_vec, func_deps, 2, 0)
  priors = rep(list(prior_dist(name="uniform",params=c(0, 1), bounds=c(0,2))),2)
  generate(mod,priors, "test_gen.stan")
  out = readChar("test_gen.stan", file.info("test_gen.stan")$size)
  expect_equal(out, REFERENCE_MODEL_1)
  
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
  mod = bp_model(e_mat, p_vec, func_deps, 5, 1)
  priors = rep(list(prior_dist(name="normal",params=c(0,.25), bounds=c(0,5))),5)
  priors[[3]] = prior_dist(name="uniform",params=c(5,15), bounds=c(5,15))
  
  generate(mod, priors, "test_gen.stan")
  out = readChar("test_gen.stan", file.info("test_gen.stan")$size)
  expect_equal(out, REFERENCE_MODEL_2)
  
  e_mat =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
  p_vec = c(1, 1, 2, 2, 1, 2)
  func_deps = c('c[1]','c[2]','c[3]', 'c[4]', 'c[5]', 'c[6]')
  mod = bp_model(e_mat, p_vec, func_deps, 6, 0)
  priors = rep(list(prior_dist(name="uniform",params=c(0,2), bounds=c(0,5))),6)
  
  generate(mod,priors, "test_gen.stan")
  out = readChar("test_gen.stan", file.info("test_gen.stan")$size)
  expect_equal(out, REFERENCE_MODEL_3)
  
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5] - c[6]/(1 + exp(c[7]*(x[1] - c[8])))') #logistic function
  mod = bp_model(e_mat, p_vec, func_deps, 8, 1)
  priors = rep(list(prior_dist(name="normal",params=c(0,.25), bounds=c(0,5))),8)
  priors[[3]] = prior_dist(name="uniform",params=c(5,15), bounds=c(5,15))
  priors[[7]] = prior_dist(name="uniform",params=c(5,15), bounds=c(5,15))
  
  generate(mod, priors, "test_gen.stan")
  out = readChar("test_gen.stan", file.info("test_gen.stan")$size)
  expect_equal(out, REFERENCE_MODEL_4)
  
  e_mat =  rbind(c(2,0,0),c(1,1,0),c(1,0,1),c(0,2,0),c(0,0,2), c(0,0,0), c(0,0,0), c(0,0,0))
  p_vec = c(1,1,1,2,3,1,2,3)
  func_deps = c('c[1]','c[2]','c[3]','c[4]','c[5]','c[6]','c[7]','c[8]')
  mod = bp_model(e_mat, p_vec, func_deps, 8, 0)
  priors = rep(list(prior_dist(name="uniform",params=c(0,2), bounds=c(0,5))),8)
  
  generate(mod, priors, "test_gen.stan")
  out = readChar("test_gen.stan", file.info("test_gen.stan")$size)
  expect_equal(out, REFERENCE_MODEL_5)
  
  file.remove("test_gen.stan")
  
  e_mat <-  matrix(c(2,0),ncol=1)
  p_vec <- c(1, 1)
  func_deps <- c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
  mod <- bp_model(e_mat, p_vec, func_deps, 5, 1)
  priors = list()
  priors[[1]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))
  priors[[2]] <- prior_dist(name="normal", params = c(.25, .25), bounds = c(0,5))
  priors[[3]] <- prior_dist(name="normal",params=c(0,2), bounds=c(0,10))
  priors[[4]] <- prior_dist(name="normal",params=c(0,2), bounds=c(-5,5))
  priors[[5]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))
  out = generate(mod, priors, simple_bd = T)
  expect_equal(out, REFERENCE_MODEL_6)
  
})