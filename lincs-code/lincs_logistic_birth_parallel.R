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
bin_means = dplyr::summarise(dplyr::group_by(small_data, drug_conc), Mean = mean(final_pop))
gr_drop = (bin_means$Mean - dplyr::lag(bin_means$Mean))[-1]
drop_idx = which.min(gr_drop)+1
gr_midpoint = bin_means[drop_idx,]$drug_conc #empirical prior

gr_range = max(bin_means[1,]$Mean - bin_means[nrow(bin_means),]$Mean,0.01) #empirical prior

ggplot() + geom_point(data=small_data, aes(x = drug_conc, y = log(final_pop/init_pop)/times))

stan_dat = create_stan_data(mod, final_pop = matrix(small_data$final_pop), init_pop = matrix(small_data$init_pop), c_mat = matrix(small_data$drug_conc), times = small_data$times)
priors = list()
priors[[1]] <- prior_dist(name="normal", params = c(0, .5), bounds = c(0,5))
priors[[2]] <- prior_dist(name="normal", params = c(gr_range, .5), bounds = c(0,5))
priors[[3]] <- prior_dist(name="normal",params=c(0,7), bounds=c(0,10))
priors[[4]] <- prior_dist(name="normal",params=c(gr_midpoint,2), bounds=c(-6,6))
priors[[5]] <- prior_dist(name="normal", params = c(0, .5), bounds = c(0,5))


generate(mod, priors, "lincs_birth_logistic.stan")

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,5), mod$nparams),ncol=2,byrow = T)
ranges[4,] <- c(gr_midpoint-1,gr_midpoint+1)
ranges[3,] <- c(0,10)
init <- uniform_initialize(ranges, 4)

if(file.exists("/michorlab/jroney/bpinference/lincs-data/compiled/lincs_birth_logistic.RDS")){
  stan_mod <- readRDS("/michorlab/jroney/bpinference/lincs-data/compiled/lincs_birth_logistic.RDS")
} else{
  stan_mod <- rstan::stan_model(file = "lincs_birth_logistic.stan")
  saveRDS(stan_mod, "/michorlab/jroney/bpinference/lincs-data/compiled/lincs_birth_logistic.RDS")
}

fit_data <- rstan::sampling(stan_mod, data = stan_dat, control = list(adapt_delta = 0.95, max_treedepth = 20), chains = 4, refresh = 1, init =init, iter=3000)

samples = data.frame(extract(fit_data))

compute_growth_curve = 
  function(params, doses){sapply(doses, function(x){params[1] + params[2]/(1 + exp(params[3]*(x - params[4]))) - params[5]})}

sample_growth_curves = apply(samples[1:800,], 1, compute_growth_curve, seq(min(small_data$drug_conc),max(small_data$drug_conc),length.out = 1000)) 
gcurves = cbind(reshape::melt(data.frame(sample_growth_curves)),dose = seq(min(small_data$drug_conc),max(small_data$drug_conc),length.out = 1000))

save(fit_data, small_data, file=sprintf("lincs-data/inference/%s_%s.rda",cell_name, drug))


warns = sprintf("%s_%s Div: %d Treedepth: %s Rhat: %s\n", cell_name, drug, check_div(fit_data), check_exceeded_treedepth(fit_data, max_depth = 20), check_rhat(fit_data))
cat(warns, file="lincs-data/inference/warnings.txt",append = TRUE)
