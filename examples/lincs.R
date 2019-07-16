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

lincs_to_stan <- function(small_data, model){ 
  small_data = dplyr::filter(lincs_data, Cell.Name == cell, Small.Molecule.Name == drug_name)
  controls = small_data$Total.Control.Cell.Count[1:max(small_data$Replicate)];
  c_mat = log(matrix(c(rep(0.001, length(controls)), small_data$Small.Mol.Concentration),ncol=1))
  final_pop = matrix(c(controls, small_data$Total.Cell.Count.After.Treatment),ncol=1)
  times = c(rep( small_data$Timepoint[1], length(controls)), small_data$Timepoint)/24
  init_pop = matrix(c(rep(small_data$Total.Cell.Count.Before.Treatment[1], length(controls)), small_data$Total.Cell.Count.Before.Treatment),ncol=1)
  return(create_stan_data(model = model, c_mat = c_mat, final_pop = final_pop, init_pop = init_pop, times= times))
}

small_data = dplyr::filter(lincs_data, Cell.Name == "MCF7", Small.Molecule.Name == "Seliciclib")
stan_dat = lincs_to_stan(small_data, mod)
priors <- rep(list(list(name="normal",params=c(0,.25), bounds=c(0,5))),5)
priors[[4]] <- list(name="uniform",params=c(-2,2), bounds=c(-2,2))
priors[[3]] <- list(name="uniform",params=c(5,15), bounds=c(5,15))

generate(mod, priors, "lincs_model.stan")

options(mc.cores = parallel::detectCores())

ranges <- matrix(rep(c(0,1), mod$nparams),ncol=2,byrow = T)
ranges[4,] <- c(-1,1)
ranges[3,] <- c(5,15)
init <- uniform_initialize(ranges, 4)

stan_mod <- rstan::stan_model(file = "lincs_model.stan")
fit_data <- rstan::sampling(stan_mod, data = stan_dat, control = list(adapt_delta = 0.95), chains = 4, refresh = 1, init =init, iter=3000)

samples = extract(fit_data,permute = FALSE)

compute_birth_curve = 
  function(params, doses){sapply(doses, function(x){params[1] + params[2]/(1 + exp(params[3]*(x - params[4])))})}

mean_birth_curve1 = mean(apply(samples[,1,1:4], 1, compute_birth_curve, seq(-5,5,length.out = 1000)))
mean_pop_curve1 = 
sample_birth_curves2 = apply(samples[1:400,2,1:4], 1, compute_birth_curve, seq(-5,5,length.out = 1000)) 
bcurves1 = cbind(reshape::melt(data.frame(sample_birth_curves1)),dose = seq(-5,5,length.out = 1000))
bcurves2 = cbind(reshape::melt(data.frame(sample_birth_curves2)),dose = seq(-5,5,length.out = 1000))

ggplot(bcurves, aes(x=dose, y=value, group = factor(variable))) + geom_line(data = bcurves1,aes(color = "chain 1"), alpha = .03) + geom_line(data = bcurves2,aes(color = "chain 2"), alpha = .03)  
