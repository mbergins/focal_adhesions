rm(list = ls())
source("/Users/mbergins/Documents/Projects/focal_adhesions/src/analyze_cell_features/FA_analysis_lib.R")
library(boot)

################################################################################
#Basic Testing Code
################################################################################
basic_test_set = c()
basic_test_set = c(NaN, exp(seq(0.1,1,length=10)), rep(1,length=10), exp(seq(1,0.1,length=10)), NaN)
basic_test_set = rbind(basic_test_set, c(NaN, NaN, exp(seq(0.1,1,length=10)), rep(1,length=10), exp(seq(1,0.1,length=10))))
basic_test_set = rbind(basic_test_set, c(exp(seq(0.1,1,length=10)), rep(1,length=10), exp(seq(1,0.1,length=10)),NaN, NaN))
basic_test_set = rbind(basic_test_set, c(rep(NaN, dim(basic_test_set)[[2]] - 1),1))

basic_test_props = data.frame(split_birth_status = rep(0,dim(basic_test_set)[[1]]), 
					    	  death_status = rep(1,dim(basic_test_set)[[1]])
				 	   		 )

basic_test_results <- gather_bilinear_models(basic_test_set,basic_test_props)

stopifnot(dim(basic_test_results$disassembly)[[1]] == dim(basic_test_set)[[1]])
stopifnot(dim(basic_test_results$assembly)[[1]] == dim(basic_test_set)[[1]])
stopifnot(all(is.na(basic_test_results$assembly$R_sq[3:4])))
stopifnot(all(is.na(basic_test_results$disassembly$R_sq[2]) & is.na(basic_test_results$disassembly$R_sq[4])))

################################################################################
#Realistic Testing Code
################################################################################
test_data <- as.matrix(read.table('../../results/focal_adhesions/time_series_01/adhesion_props/lin_time_series/Average_adhesion_signal.csv', header=FALSE, sep=','));
filt_data <- test_data[is.nan(test_data[,1]) & is.nan(test_data[, dim(test_data)[[2]] ]), ]

min_max = c()
for (i in 1:dim(filt_data)[[1]]) {
	line = na.omit(filt_data[i,])
	
	if (is.nan(filt_data[i,1]) & is.nan(filt_data[i,dim(filt_data)[[2]]]) & length(line) > 20) {
		min_max = rbind(min_max, c(min(line), max(line)))
	}
}
min_max = data.frame(low = min_max[,1], high = min_max[,2])

as_time = 15;
stable_lifetime = 20;
error_amounts = seq(0.001,0.008,length=10)
data_points = list()
for (i in error_amounts[1:4]) {
	data_sets = c()
	for (j in 1:200) {
		h_l = min_max[sample(1:dim(min_max)[[1]],1),]
		
		
		k = log(h_l$high/h_l$low)/as_time;
		assembly_phase = h_l$low*exp(k * seq(0, as_time - 1)) + rnorm(as_time,sd=i);
		
		data = c(NaN, 
				 assembly_phase,
				 rep(h_l$high,stable_lifetime) + rnorm(stable_lifetime,sd=0.01),
				 exp(seq(h_l$high,h_l$low,length=10))*(h_l$high/exp(h_l$high)), 
				 NaN)
-		data_sets = rbind(data_sets, data);
	}
	data_props = data.frame(split_birth_status = rep(0,dim(data_sets)[[1]]), 
						    death_status = rep(1,dim(data_sets)[[1]])
					 	   )
	these_results = gather_bilinear_models(data_sets, data_props,debug=1)
#	print(shapiro.test(these_results$assembly$length))
#	hist(these_results$assembly$length)
	boot_samp = boot(these_results$assembly$length, function(data,indexes) mean(data[indexes], na.rm=T), 10000)
	boot_conf = boot.ci(boot_samp,type="perc")
	data_points$upper = c(data_points$upper, boot_conf$perc[5])
	data_points$lower = c(data_points$lower, boot_conf$perc[4])
	data_points$mean = c(data_points$mean, boot_samp$t0)

	boot_samp = boot(these_results$assembly$R_sq, function(data,indexes) mean(data[indexes], na.rm=T),10000)
	boot_conf = boot.ci(boot_samp,type="perc")
	data_points$R_upper = c(data_points$R_upper, boot_conf$perc[5])
	data_points$R_lower = c(data_points$R_lower, boot_conf$perc[4])
	data_points$R_mean = c(data_points$R_mean, boot_samp$t0)
	print(paste('Error Level:', i, ' done'))
}

library(Hmisc)