context("test that full stan models are properly generated")

test_that("incorrect models are not generated", {
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  func_deps = c('c[1]','c[2]')
  mod = bp_model(e_mat, p_vec, func_deps, 2, 0)
  priors = rep(list(prior_dist(name="uniform",params=c(0, 1), bounds=c(0,2))),1)
  expect_error(generate(mod,priors), "incorrect number of priors for model!")
  
  priors = rep(list(prior_dist(name="uniform",params=c(0, 1), bounds=c(0,2))),3)
  expect_error(generate(mod,priors), "incorrect number of priors for model!")
})
  
test_that("generated models match known reference models", {
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  func_deps = c('c[1]','c[2]')
  mod = bp_model(e_mat, p_vec, func_deps, 2, 0)
  priors = rep(list(prior_dist(name="uniform",params=c(0, 1), bounds=c(0,2))),2)
  code = generate(mod,priors)
  expect_equal(digest::digest(code, "md5"), "326cc93ef56ade2de45324dea009237f")
  
  e_mat =  matrix(c(2,0),ncol=1)
  p_vec = c(1, 1)
  func_deps = c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
  mod = bp_model(e_mat, p_vec, func_deps, 5, 1)
  priors = rep(list(prior_dist(name="normal",params=c(0,.25), bounds=c(0,5))),5)
  priors[[3]] = prior_dist(name="uniform",params=c(5,15), bounds=c(5,15))
  
  code = generate(mod, priors)
  expect_equal(digest::digest(code, "md5"), "d91e460da4f4d43dc3896d93a330b29f")
  
  e_mat =  rbind(c(2, 0),c(1,1),c(1,1), c(0,2), c(0,0),c(0,0))
  p_vec = c(1, 1, 2, 2, 1, 2)
  func_deps = c('c[1]','c[2]','c[3]', 'c[4]', 'c[5]', 'c[6]')
  mod = bp_model(e_mat, p_vec, func_deps, 6, 0)
  priors = rep(list(prior_dist(name="uniform",params=c(0,2), bounds=c(0,5))),6)
  
  code = generate(mod,priors)
  expect_equal(digest::digest(code, "md5"), "07da76feb881304b67ff4e2efc9531a8")
  
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
  code = generate(mod, priors, simple_bd = T)
  expect_equal(digest::digest(code, "md5"), "273308f4dfd7b2e7a0bd06e263a5d7da")
  
  e_mat <-  rbind(c(2,0),c(0,1),c(0,0))
  p_vec <- c(1, 1, 2)

  func_deps <- c('c[1]','c[2]', 'c[3]')
  priors <- rep(list(list(name="normal",params=c(0, .25), bounds=c(0,2))),3)
  
  mod <- bp_model(e_mat, p_vec, func_deps, 3, 0)
  code <- generate(mod, priors, noise_model = T)
  expect_equal(digest::digest(code, "md5"), "feca3712f0ff8821011f46942151cbd5")
})