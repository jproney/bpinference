# run this on the DFCI Kraken cluster
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all()
cell_name = args[1]
drug = args[2]

delete_cols = c("Timepoint.Unit","Small.Mol.Conc.Unit","LJP.Library.Plate","Assay.Well","Normalized.Growth.Rate.Inhibition.Value","Cell.HMS.LINCS.ID","Small.Molecule.HMS.LINCS.ID","Relative.Cell.Count")
lincs_data = read.csv("lincs-data/dataset_20245_20181218181013.csv")
lincs_data = lincs_data[,!names(lincs_data)%in%delete_cols]

e_mat <-  matrix(c(2,0),ncol=1)
p_vec <- c(1, 1)
func_deps <- c('c[1] + c[2]/(1 + exp(c[3]*(x[1] - c[4])))','c[5]') #logistic function
mod <- bp_model(e_mat, p_vec, func_deps, 5, 1)

prepare_data <- function(cellname, drug){ 
  small_data = dplyr::filter(lincs_data, Cell.Name == cellname, Small.Molecule.Name == drug)
  controls = small_data$Total.Control.Cell.Count[1:max(small_data$Replicate)];
  c_mat = log(c(rep(0.001, length(controls)), small_data$Small.Mol.Concentration))
  final_pop = c(controls, small_data$Total.Cell.Count.After.Treatment)
  times = c(rep( small_data$Timepoint[1], length(controls)), small_data$Timepoint)/24
  init_pop = c(rep(small_data$Total.Cell.Count.Before.Treatment[1], length(controls)), small_data$Total.Cell.Count.Before.Treatment)
  new_data = data.frame(cbind(drug_conc = c_mat, final_pop = final_pop, times = times, init_pop = init_pop))
  return(new_data)
}

small_data = prepare_data(cell_name, drug)
bin_means = dplyr::summarise(dplyr::group_by(small_data, drug_conc), Mean = mean(log(final_pop/init_pop)/3))
gr_drop = (bin_means$Mean - dplyr::lag(bin_means$Mean))[-1]
drop_idx = which.min(gr_drop)+1
gr_midpoint = bin_means[drop_idx,]$drug_conc #empirical prior

gr_range = max(bin_means[1,]$Mean - bin_means[nrow(bin_means),]$Mean,0.01) #empirical prior

#Theta1 ~ normal(0,0.25);
#Theta2 ~ normal(0.629613491882936,0.5);
#Theta3 ~ uniform(1,10);
#Theta4 ~ uniform(-3,3);
#Theta5 ~ normal(0,0.25);

#real<lower=0, upper=5> Theta1;
#real<lower=0, upper=5> Theta2;
#real<lower=1, upper=10> Theta3;
#real<lower=-3, upper=3> Theta4;
#real<lower=0, upper=5> Theta5;

stan_dat = create_stan_data(mod, final_pop = matrix(small_data$final_pop), init_pop = matrix(small_data$init_pop), c_mat = matrix(small_data$drug_conc), times = small_data$times)
priors = list()
priors[[1]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))
priors[[2]] <- prior_dist(name="normal", params = c(gr_range, .5), bounds = c(0,5))
priors[[3]] <- prior_dist(name="uniform",params=c(1,10), bounds=c(1,10))
priors[[4]] <- prior_dist(name="normal",params=c(gr_midpoint,1), bounds=c(-6,6))
priors[[5]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))


model_str = generate(mod, priors) #regenerate every time

options(mc.cores = parallel::detectCores())

stan_mod <- rstan::stan_model(model_code = model_str)

fit_data <- rstan::sampling(stan_mod, data = stan_dat, control = list(adapt_delta = 0.95, max_treedepth = 20), chains = 4, refresh = 1, iter=3000)

samples = data.frame(extract(fit_data))

save(fit_data, small_data, file=sprintf("lincs-data/inference/%s_%s.rda",cell_name, drug))


warns = sprintf("%s_%s Div: %d Treedepth: %s Rhat: %s\n", cell_name, drug, check_div(fit_data), check_exceeded_treedepth(fit_data, max_depth = 20), check_rhat(fit_data))
cat(warns, file="lincs-data/inference/warnings.txt",append = TRUE)
