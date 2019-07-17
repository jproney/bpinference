delete_cols = c("Timepoint.Unit","Small.Mol.Conc.Unit","LJP.Library.Plate","Assay.Well","Normalized.Growth.Rate.Inhibition.Value","Cell.HMS.LINCS.ID","Small.Molecule.HMS.LINCS.ID","Relative.Cell.Count")
lincs_data = read.csv("lincs-data/dataset_20245_20181218181013.csv")
lincs_data = lincs_data[,!names(lincs_data)%in%delete_cols]

cell_lines = c("MCF7","SK-BR-3","MDA-MB-231 BT-20","Hs 578T","MCF 10A")
drugs =   c('Seliciclib','AT-7519','AZD7762','AZD8055','Sorafenib','CP466722','CP724714','Alvocidib','GW843682X','HG-5-113-01','HG-5-88-01','HG-6-64-01','Neratinib','JW-7-24-1',
            'Dasatinib','Tozasertib','Imatinib','NVP-TAE684','CGP60474','PD173074','Crizotinib','BMS345541','KIN001-043','Saracatinib','Sigma A6730','WH-4-025','R406','BI-2536',
            'A443654','SB590885','Pictilisib','PD184352','PLX-4720','Lapatinib','Sirolimus','ZSTK474','BX-912','Selumetinib','MK2206','AZD-6482','NU7441','OSI-027','WYE-125132',
            'Barasertib','Vemurafenib','Enzastaurin','Palbociclib','PF562271','PHA-793887','QL-X-138','QL-XII-47','Torin1','Torin2','WZ-4-145','WZ3105','XMD11-50','XMD11-85h',
            'XMD16-144','Erlotinib','Gefitinib','Nilotinib','JNK-9L','PD0325901','Geldanamycin','YM 201636','TWS119','PF477736','LDN-193189','PF431396','Celastrol','SU11274',
            'Canertinib','NVP-AEW541','PHA-665752','PI103','Dovitinib','GSK 690693','SNS-032','Afatinib','GSK1904529A','Linsitinib','Ruxolitinib','Momelotinib','Fedratinib','Trametinib',
            'Omipalisib','Buparlisib','XL147','Y39983','Nintedanib','Foretinib','AZD 5438','Pelitinib','Luminespib','AZD8330','TGX221','GSK1059615','Brivanib','CHIR-99021','Linifanib',
            'PIK-93','Dactolisib','Alpelisib','GDC-0980','Mitoxantrone','Radicicol','Withaferin A')


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

cell_name = "Hs-578T"
drug = 'YM-201636'
small_data = prepare_data(cell_name, drug)
gr = log(small_data$final_pop/small_data$init_pop)/3
gr_range = mean(gr[1:3]) - mean(tail(gr,3)) #empirical prior mean

ggplot() + geom_point(data=small_data, aes(x = drug_conc, y = log(final_pop/init_pop)/times))

stan_dat = create_stan_data(mod, final_pop = matrix(small_data$final_pop), init_pop = matrix(small_data$init_pop), c_mat = matrix(small_data$drug_conc), times = small_data$times)
priors = list()
priors[[1]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))
priors[[2]] <- prior_dist(name="normal", params = c(gr_range, .25), bounds = c(0,5))
priors[[3]] <- prior_dist(name="uniform",params=c(1,10), bounds=c(1,10))
priors[[4]] <- prior_dist(name="uniform",params=c(-3,3), bounds=c(-3,3))
priors[[5]] <- prior_dist(name="normal", params = c(0, .25), bounds = c(0,5))

generate(mod, priors, "lincs_birth_logistic.stan")

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1), mod$nparams),ncol=2,byrow = T)
ranges[4,] <- c(-1,1)
ranges[3,] <- c(1,10)
init <- uniform_initialize(ranges, 4)

stan_mod <- rstan::stan_model(file = "lincs_birth_logistic.stan")
fit_data <- rstan::sampling(stan_mod, data = stan_dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init, iter=3000)

samples = data.frame(extract(fit_data))

compute_growth_curve = 
  function(params, doses){sapply(doses, function(x){params[1] + params[2]/(1 + exp(params[3]*(x - params[4]))) - params[5]})}

sample_growth_curves = apply(samples[1:800,], 1, compute_growth_curve, seq(min(small_data$drug_conc),max(small_data$drug_conc),length.out = 1000)) 
gcurves = cbind(reshape::melt(data.frame(sample_growth_curves)),dose = seq(min(small_data$drug_conc),max(small_data$drug_conc),length.out = 1000))

save(samples, small_data, file=sprintf("lincs-data/inference/%s_%s.rda",cell_name, drug))
png(sprintf("lincs-data/inference/%s_%s.png",cell_name, drug))
plt <- ggplot() + geom_line(data=gcurves, aes(x=dose, y=value, group = factor(variable)),alpha=.03, color="red") + geom_point(data=small_data, aes(x = drug_conc, y = log(final_pop/init_pop)/times))
print(plt)
dev.off()


warns = sprintf("%s_%s Div: %d Treedepth: %s Rhat: %s\n", cell_name, drug, check_div(fit_data), check_exceeded_treedepth(fit_data), check_rhat(fit_data))
cat(warns, file="lincs-data/inference/warnings.txt",append = TRUE)
