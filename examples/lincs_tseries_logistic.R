lincs = read.csv("lincs-data/timeseries/lincs_ts.csv")
lincs <- lincs %>% group_by(AssayWell, Plate, Concentration) %>% mutate(PrevCellCount = dplyr::lag(CellCountAfterTreatment), DeltaT = Timepoint - dplyr::lag(Timepoint))
lincs <- lincs %>% filter(!is.na(PrevCellCount))
drug_data <- lincs %>% filter((is.na(Drug) | Drug == "Alpelisib") & Timepoint < 40)

func_deps <- c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
ndep = 1
nparam = 5
mod = bp_model_simple_birth_death(func_deps, nparam, ndep)

priors[[1]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))
priors[[2]] <- prior_dist(name="normal", params = c(0, .5), bounds = c(0,5))
priors[[3]] <- prior_dist(name="uniform",params=c(1,10), bounds=c(1,10))
priors[[4]] <- prior_dist(name="normal",params=c(0,1), bounds=c(-6,6))
priors[[5]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))

generate(mod, priors, "lincs_ts.stan", simple_bd = T)

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nrow(e_mat)),nrow(e_mat),2,byrow = T)
init <- uniform_initialize(ranges, 4)

drug_data[drug_data$Concentration == 0,]$Concentration = .001
dat <- create_stan_data(mod, matrix(drug_data$CellCountAfterTreatment, ncol=1), matrix(drug_data$PrevCellCount,ncol=1), drug_data$DeltaT, c_mat = matrix(log(drug_data$Concentration),ncol=1), simple_bd = T)

stan_mod <- rstan::stan_model(file = "lincs_ts.stan")
fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, iter = 3000, warmup = 1000)
file.remove("lincs_ts.stan")