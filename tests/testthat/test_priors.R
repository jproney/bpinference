context("test instantiation of prior objects")

test_that("Invalid prior objects are not instantiated",{
  expect_error(prior_dist(name = "are you", params = c(1,2,3), bounds = c(2,3)), "Invalid prior name!")
  expect_error(prior_dist(name = 24, params = c(1,2,3), bounds = c(2,3)), "prior name must be a string!")
  expect_error(prior_dist(name = "normall", params = c(1,2,3), bounds = c(2,3)), "Invalid prior name!")
  expect_error(prior_dist(name = "normal", params = c(1,2,3), bounds = c(2,3)), "Incorrect number of parameters for prior distribution")
  expect_error(prior_dist(name = "normal", params = c(1,0), bounds = c(2,3)), "Positive parameter required here!")
  expect_error(prior_dist(name = "normal", params = c(1,-3), bounds = c(2,3)), "Positive parameter required here!")
  expect_error(prior_dist(name = "lognormal", params = c(1,-3), bounds = c(2,3)), "Positive parameter required here!")
  expect_error(prior_dist(name = "exponential", params = c(1,-3), bounds = c(2,3)), "Incorrect number of parameters for prior distribution")
  expect_error(prior_dist(name = "exponential", params = c(0), bounds = c(2,3)), "Positive parameter required here!")
  expect_error(prior_dist(name = "exponential", params = c(-1), bounds = c(2,3)), "Positive parameter required here!")
  expect_error(prior_dist(name = "exponential", params = "feeling", bounds = c(2,3)), "parameters should be numeric!")
  expect_error(prior_dist(name = "gamma", params = c(1,1), bounds = "it"), "bounds should be a list of 2 numbers!")
  expect_error(prior_dist(name = "gamma", params = c(1,1), bounds = "now"), "bounds should be a list of 2 numbers!")
  expect_error(prior_dist(name = "gamma", params = c(1,1), bounds = c(4,5,6)), "bounds should be a list of 2 numbers!")
  expect_error(prior_dist(name = "gamma", params = c(1,1), bounds = c(25,24)), "Error: lower bounder greater or equal to upper bound!")
  expect_error(prior_dist(name = "gamma", params = c(-1,1), bounds = c(24,25)), "Positive parameter required here!")
  expect_error(prior_dist(name = "uniform", params = c(1,-1), bounds = c(-2,2)), "upper bound must be greater than lower in uniform distribution!")
  expect_error(prior_dist(name = "MR KRABS?!?!?!?", params = c(1,2,3), bounds = c(2,3)), "Invalid prior name!")
})

test_that("Correct priors are parsed properly", {
  expect_equal(parse_prior(prior_dist(name = "normal", params = c(0,1), bounds = c(-1,1)),1)[1], "\tTheta1 ~ normal(0,1);\n")
  expect_equal(parse_prior(prior_dist(name = "normal", params = c(0,1), bounds = c(-1,1)),1)[2], "\treal<lower=-1, upper=1> Theta1;\n")
  expect_equal(parse_prior(prior_dist(name = "uniform", params = c(0,2), bounds = c(0,3)),4)[1], "\tTheta4 ~ uniform(0,2);\n")
  expect_equal(parse_prior(prior_dist(name = "uniform", params = c(0,2), bounds = c(0,3)),4)[2], "\treal<lower=0, upper=3> Theta4;\n")
  expect_equal(parse_prior(prior_dist(name = "normal", params = c(0,.4324)),1)[1], "\tTheta1 ~ normal(0,0.4324);\n")
  expect_equal(parse_prior(prior_dist(name = "normal", params = c(0,.4324)),1)[2], "\treal Theta1;\n")
  expect_equal(parse_prior(prior_dist(name = "gamma", params = c(.5,5.77), bounds = c(0,58.18)),1)[1], "\tTheta1 ~ gamma(0.5,5.77);\n")
  expect_equal(parse_prior(prior_dist(name = "gamma", params = c(.5,5.77), bounds = c(0,58.18)),1)[2], "\treal<lower=0, upper=58.18> Theta1;\n")
})