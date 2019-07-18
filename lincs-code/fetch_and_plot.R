cell_lines = c("MCF7","SK-BR-3","MDA-MB-231 BT-20","Hs 578T","MCF 10A")
drugs =   c('Seliciclib','AT-7519','AZD7762','AZD8055','Sorafenib','CP466722','CP724714','Alvocidib','GW843682X','HG-5-113-01','HG-5-88-01','HG-6-64-01','Neratinib','JW-7-24-1',
            'Dasatinib','Tozasertib','Imatinib','NVP-TAE684','CGP60474','PD173074','Crizotinib','BMS345541','KIN001-043','Saracatinib','Sigma A6730','WH-4-025','R406','BI-2536',
            'A443654','SB590885','Pictilisib','PD184352','PLX-4720','Lapatinib','Sirolimus','ZSTK474','BX-912','Selumetinib','MK2206','AZD-6482','NU7441','OSI-027','WYE-125132',
            'Barasertib','Vemurafenib','Enzastaurin','Palbociclib','PF562271','PHA-793887','QL-X-138','QL-XII-47','Torin1','Torin2','WZ-4-145','WZ3105','XMD11-50','XMD11-85h',
            'XMD16-144','Erlotinib','Gefitinib','Nilotinib','JNK-9L','PD0325901','Geldanamycin','YM 201636','TWS119','PF477736','LDN-193189','PF431396','Celastrol','SU11274',
            'Canertinib','NVP-AEW541','PHA-665752','PI103','Dovitinib','GSK 690693','SNS-032','Afatinib','GSK1904529A','Linsitinib','Ruxolitinib','Momelotinib','Fedratinib','Trametinib',
            'Omipalisib','Buparlisib','XL147','Y39983','Nintedanib','Foretinib','AZD 5438','Pelitinib','Luminespib','AZD8330','TGX221','GSK1059615','Brivanib','CHIR-99021','Linifanib',
            'PIK-93','Dactolisib','Alpelisib','GDC-0980','Mitoxantrone','Radicicol','Withaferin A')

cell_line = "SK-BR-3"
drug = "Torin1"

command = sprintf("scp jamesr@kraken.dfci.harvard.edu:/michorlab/jroney/bpinference/lincs-data/inference/%s_%s.rda %s_%s.rda", cell_line, drug, cell_line, drug)
system(command, wait = TRUE)

load(sprintf("%s_%s.rda", cell_line, drug))
#samples = data.frame(extract(fit_data))
samples = extract(fit_data, permute=F)
compute_growth_curve = 
  function(params, doses){sapply(doses, function(x){params[1] + params[2]/(1 + exp(params[3]*(x - params[4]))) - params[5]})}

sample_growth_curves = apply(samples[1:800,1,], 1, compute_growth_curve, seq(min(small_data$drug_conc),max(small_data$drug_conc),length.out = 1000)) 
gcurves = cbind(reshape::melt(data.frame(sample_growth_curves)),dose = seq(min(small_data$drug_conc),max(small_data$drug_conc),length.out = 1000))

ggplot() + geom_line(data=gcurves, aes(x=dose, y=value, group = factor(variable)),alpha=.03, color="red") + geom_point(data=small_data, aes(x = drug_conc, y = log(final_pop/init_pop)/times))

bin_means = dplyr::summarise(dplyr::group_by(small_data, drug_conc), Mean = mean(final_pop))
gr_drop = (bin_means$Mean - dplyr::lag(bin_means$Mean))[-1]
drop_idx = which.min(gr_drop)+1
heuristic_midpoint = bin_means[drop_idx,]$drug_conc