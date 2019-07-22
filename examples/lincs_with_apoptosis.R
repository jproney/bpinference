lincs_data = read.csv("lincs-data/dataset_20217.csv")
cont = dplyr::filter(lincs_data, Drug.Concentration==0, Cell.Line.Name == "C32")

init_live = c(dplyr::filter(cont,Time.Point == 24)$Live.Cells, dplyr::filter(cont,Time.Point == 48)$Live.Cells)
init_dead = c(dplyr::filter(cont,Time.Point == 24)$Number.of.Apoptotic.Cells, dplyr::filter(cont,Time.Point == 48)$Number.of.Apoptotic.Cells)

final_live = c(dplyr::filter(cont,Time.Point == 48)$Live.Cells, dplyr::filter(cont,Time.Point == 72)$Live.Cells)
final_dead = c(dplyr::filter(cont,Time.Point == 48)$Number.of.Apoptotic.Cells, dplyr::filter(cont,Time.Point == 72)$Number.of.Apoptotic.Cells)

init_pop = cbind(init_live, init_dead)
final_pop = cbind(final_live, final_dead)
times = rep(1, nrow(init_pop))

e_mat = rbind(c(2,0),c(0,1), c(0,0))
p_vec = c(1,1,2)
func_deps = c("c[1]","c[2]","c[3]")
ndep = 0
nparam = 3
mod = bp_model(e_mat, p_vec, func_deps, nparam, ndep)

priors <- rep(list(list(name="normal",params=c(0, .25), bounds=c(0,5))),3)

generate(mod, priors, "lincs_apop.stan")

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init <- uniform_initialize(ranges, 4)

dat <- create_stan_data(mod, final_pop, init_pop, times)

stan_mod <- rstan::stan_model(file = "lincs_apop.stan")
fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, iter = 3000, warmup = 1000)
file.remove("one_type.stan")