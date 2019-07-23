lincs = read.csv("lincs-data/timeseries/lincs_ts.csv")
lincs <- lincs %>% group_by(AssayWell, Plate) %>% mutate(PrevCellCount = dplyr::lag(CellCountAfterTreatment), DeltaT = Timepoint - dplyr::lag(Timepoint))
lincs <- lincs %>% filter(!is.na(PrevCellCount))
controls <- lincs %>% filter(is.na(Drug) & Timepoint < 40)

ggplot(controls, aes(x = Timepoint, y = CellCountAfterTreatment)) + geom_point(aes(color = AssayWell))

func_deps = c("c[1]","c[2]")
ndep = 0
nparam = 2
mod = bp_model_simple_birth_death(func_deps, nparam, ndep)

priors <- rep(list(list(name="normal",params=c(0, .05), bounds=c(0,5))),2)

generate(mod, priors, "lincs_ts.stan", simple_bd = T)

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nparam),nparam,2,byrow = T)
init <- uniform_initialize(ranges, 4)

dat <- create_stan_data(mod, controls$CellCountAfterTreatment, controls$PrevCellCount, controls$DeltaT, simple_bd = T)

stan_mod <- rstan::stan_model(file = "lincs_ts.stan")
fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, iter = 3000, warmup = 1000)
file.remove("lincs_ts.stan")