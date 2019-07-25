lincs1 = read.csv("lincs-data/timeseries/lincs_ts.csv")
lincs2 = read.csv("lincs-data/timeseries/lincs_tseries2.csv")
lincs3 = read.csv("lincs-data/timeseries/lincs_tseries3.csv")
lincs = read.csv("lincs-data/timeseries/lincs-ts-master.csv")
lincs = lincs1
lincs <- lincs %>% group_by(AssayWell, Plate) %>% mutate(PrevCellCount = dplyr::lag(CellCountAfterTreatment), DeltaT = Timepoint - dplyr::lag(Timepoint))
lincs <- lincs %>% filter(!is.na(PrevCellCount))
controls <- lincs %>% filter(is.na(Drug)  & Timepoint < 40 &!any(CellCountAfterTreatment < 500))

ggplot(controls, aes(x = Timepoint, y = log(CellCountAfterTreatment))) + geom_point(aes(color = AssayWell)) + ggtitle("Growth of Control MCF10A Cells") + ylab("Cell Count")

mu_hat_c = sum(controls$CellCountAfterTreatment)/sum(controls$PrevCellCount)
sig_hat_c = sum((controls$CellCountAfterTreatment - mu_hat_c*controls$PrevCellCount)**2)/sum(controls$PrevCellCount)
norm = dnorm(seq(-3.5,3.5,.01), 0, 1)
norm = norm/max(norm)*90
dens = data.frame(cbind(norm = norm, x = seq(-3.5,3.5,.01)))

controls <- controls %>% mutate(resid = (CellCountAfterTreatment - mu_hat_c*PrevCellCount)/(sqrt(sig_hat_c*PrevCellCount)))
ggplot(controls , aes(x = resid)) + geom_histogram(binwidth = .2,fill = "white", color = "black") + geom_line(data = dens, aes(x=x, y = norm)) + xlab("Rescaled Residual Error From Weighted Least Squares Regression") + ylab("Number of Occurances") + ggtitle("Approximate Normality of in vitro Cell Growth")

drug <- lincs %>% filter(Drug == "Dasatinib" & !any(CellCountAfterTreatment < 500))
conc = unique(drug$Concentration)
limits = matrix(c(0,50,15,50,0,50,0,40,0,50,10,50,10,50,15,50,20,60,20,60,20,60,20,60),ncol = 2, byrow = T)
mu_hat = rep(0,12)
sig_hat = rep(0,12)

i = 1
for(i in 1:12){
drug_at_conc = drug %>% filter(Concentration == conc[i] & limits[i,1] < Timepoint & Timepoint < limits[i,2])
mu_hat[i] = sum(drug_at_conc$CellCountAfterTreatment)/sum(drug_at_conc$PrevCellCount)
sig_hat[i] = sum((drug_at_conc$CellCountAfterTreatment - mu_hat[i]*drug_at_conc$PrevCellCount)**2)/sum(drug_at_conc$PrevCellCount)
ggplot(drug_at_conc, aes(x = Timepoint, y = (CellCountAfterTreatment))) + geom_point(aes(color = AssayWell)) + ylab("Cell Count") + ggtitle(sprintf("Dasatinib %f", conc[i]))
}

drug <- lincs %>% filter(Drug == "Etoposide" & !any(CellCountAfterTreatment < 500))
conc = unique(drug$Concentration)
limits = matrix(c(0,50,10,45,0,40,0,40,10,50,0,50,10,50,0,30,0,30),ncol = 2, byrow = T)
mu_hat = rep(0,9)
sig_hat = rep(0,9)

for(i in 1:9){
drug_at_conc = drug %>% filter(Concentration == conc[i] & limits[i,1] < Timepoint & Timepoint < limits[i,2])
mu_hat[i] = sum(drug_at_conc$CellCountAfterTreatment)/sum(drug_at_conc$PrevCellCount)
sig_hat[i] = sum((drug_at_conc$CellCountAfterTreatment - mu_hat[i]*drug_at_conc$PrevCellCount)**2)/sum(drug_at_conc$PrevCellCount)
ggplot(drug_at_conc, aes(x = Timepoint, y = log(CellCountAfterTreatment))) + geom_point(aes(color = AssayWell))
}

drug <- lincs %>% filter(Drug == "Taxol" & !any(CellCountAfterTreatment < 500))
conc = unique(drug$Concentration)
limits = matrix(c(0,45,0,40,10,45,0,40,10,40,0,50,10,50,0,30,0,30),ncol = 2, byrow = T)
mu_hat = rep(0,6)
sig_hat = rep(0,6)

#for(i in 1:6){
i=1
drug_at_conc = drug %>% filter(Concentration == conc[i] & limits[i,1] < Timepoint & Timepoint < limits[i,2])
mu_hat[i] = sum(drug_at_conc$CellCountAfterTreatment)/sum(drug_at_conc$PrevCellCount)
sig_hat[i] = sum((drug_at_conc$CellCountAfterTreatment - mu_hat[i]*drug_at_conc$PrevCellCount)**2)/sum(drug_at_conc$PrevCellCount)
ggplot(drug_at_conc, aes(x = Timepoint, y = log(CellCountAfterTreatment))) + geom_point(aes(color = AssayWell)) + ggtitle(conc[i])
#}
#drug <- drug %>% mutate(resid = (CellCountAfterTreatment - mu_hat[i]*PrevCellCount)/(sqrt(sig_hat[i]*PrevCellCount)))

plot(c(mu_hat_c, mu_hat))
plot(c(sig_hat_c, sig_hat))

func_deps = c("c[1]","c[2]")
ndep = 0
nparam = 2
mod = bp_model_simple_birth_death(func_deps, nparam, ndep)

priors <- rep(list(list(name="normal",params=c(0, .25), bounds=c(0,5))),2)

generate(mod, priors, "lincs_ts.stan", simple_bd = T)

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1),nparam),nparam,2,byrow = T)
init <- uniform_initialize(ranges, 4)


get_group = function(x,n) x %>% nest %>% slice(n) %>% unnest(data)
n_batch = n_groups(controls)
pop_vec = matrix(rep(0, 28*50),28,50)
init_pop = matrix(rep(0, 28*50),28,50)
times = matrix(rep(0, 28*50),28,50)
batch_sizes = rep(0, n_batch)
for(i in 1:n_groups(controls)){
  g =  get_group(controls,i)
  l = nrow(g)
  pop_vec[i,] = c(g$CellCountAfterTreatment, rep(0,50-l))
  init_pop[i,] =  c(g$PrevCellCount, rep(0,50-l))
  times[i,] =  c(g$DeltaT, rep(0,50-l))
  batch_sizes[i] = l
}

init_pop = init_pop[,1:max(batch_sizes)]
pop_vec = pop_vec[,1:max(batch_sizes)]
times = times[,1:max(batch_sizes)]

dat = list(nbatch = n_batch, max_datapts = max(batch_sizes), ndatapts = batch_sizes, pop_vec = pop_vec, init_pop = init_pop, times = times)
stan_mod <- rstan::stan_model(file = "lincs_ts_batch.stan")
fit_data <- rstan::sampling(stan_mod, data = dat, control = list(adapt_delta = 0.95), chains = 1, refresh = 1, iter = 3000, warmup = 1000)
file.remove("lincs_ts.stan")