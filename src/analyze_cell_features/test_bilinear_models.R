rm(list = ls())
source("bilinear_modeling.R")
library(boot)

################################################################################
#Basic Testing Code
################################################################################

# best_model_results = list();
# 
# assembly_phase = exp(seq(log(0.25),log(0.8),length=10));
# stability_phase = rep(0.8,length=10);
# disassembly_phase = exp(seq(log(0.8),log(0.25),length=10));
# 
# for (i in 1:1000) {
#     full_exp = c(assembly_phase + rnorm(length(assembly_phase),mean=0,sd=sd(assembly_phase)/3),
#         stability_phase + rnorm(length(stability_phase),mean=0,sd=sd(assembly_phase)/3),
#         disassembly_phase + rnorm(length(disassembly_phase),mean=0,sd=sd(disassembly_phase)/3));
#     
#     data_set = list(value = full_exp,time=seq(0,along.with=full_exp,by=1));
#     assembly_model_props = fit_all_possible_log_models(data_set);
#     assembly_model_props$type = "assembly";
# 
#     data_set = list(value = rev(full_exp),time=seq(0,along.with=full_exp,by=1));
#     disassembly_model_props = fit_all_possible_log_models(data_set);
#     disassembly_model_props$type = "disassembly";
#     
#     best_indexes = find_optimum_model_indexes(assembly_model_props,disassembly_model_props);
# 
#     best_assembly = assembly_model_props[best_indexes[1],];
#     best_disassembly = disassembly_model_props[best_indexes[2],];
#     
#     best_model_results = rbind(best_model_results,best_assembly,best_disassembly);
#     if (i %% 100 == 0) print(i)
# }

# data_set = matrix(runif(1000*100),1000,100);
# exp_props = list(split_birth_status = rep(0,dim(data_set)[1]),
#     death_status = rep(1,dim(data_set)[1]));
# no_hits_models = build_bilinear_models(data_set,exp_props);

################################################################################
#Realistic Testing Code
################################################################################

rand_data_set = matrix(NA,1000,100);

for (i in 1:dim(rand_data_set)[1]) {
    min_val = runif(1,min=0.1,max=0.5);
    max_val = runif(1,min=0.5,max=0.9);

    phase_lengths = round(runif(3)*10)+10;

    assembly_phase = exp(seq(log(min_val), log(max_val),length=phase_lengths[1]));
    stability_phase = rep(max_val,length=phase_lengths[2]);
    disassembly_phase = exp(seq(log(max_val),log(min_val),length=phase_lengths[3]));

    assembly_phase = assembly_phase + rnorm(length(assembly_phase),mean=0,sd=sd(assembly_phase)/3);
    stability_phase = stability_phase + rnorm(length(stability_phase),mean=0,sd=sd(assembly_phase)/3);
    disassembly_phase = disassembly_phase + rnorm(length(disassembly_phase),mean=0,sd=sd(disassembly_phase)/3);

    full_set = c(assembly_phase,stability_phase,disassembly_phase);
    rand_data_set[i,2:(length(full_set)+1)] = full_set;
}

exp_props = list(split_birth_status = rep(0,dim(rand_data_set)[1]),
    death_status = rep(1,dim(rand_data_set)[1]));

# sample_models = build_bilinear_models(rand_data_set,exp_props,min.phase.length=10,time.spacing=2.5)

# dir = '~/Documents/Projects/focal_adhesions/trunk/results/emma/2xKD/Pax_01/adhesion_props/';
dir = '~/Documents/Projects/focal_adhesions/trunk/results/focal_adhesions/time_series_04/adhesion_props/';

actual_data_set = as.matrix(read.csv(file.path(dir,'lin_time_series/Average_adhesion_signal.csv'),header=F));
actual_exp_props = read.csv(file.path(dir,'single_lin.csv'));

actual_models = build_bilinear_models(actual_data_set, actual_exp_props);
actual_models$exp_props = actual_exp_props

rate_filters = produce_rate_filters(actual_models);

write_assembly_disassembly_periods(actual_models,dir)
