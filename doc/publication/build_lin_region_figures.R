rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Result loading
################################################################################
raw_data <- list();

for (i in 5:9) {
	exp_dirs <- Sys.glob(file.path('../../results/lin_region*/FA/',i,'*/adhesion_props/models/'));
	exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
	raw_data$wild_type[[paste("results_",i,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));
	
	exp_dirs <- Sys.glob(file.path('../../results/lin_region*/S178A/',i,'*/adhesion_props/models/'));
	exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
	raw_data$S178A[[paste("results_",i,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));
}

exp_dirs <- Sys.glob(file.path('../../results/focal_adhesions/*/adhesion_props/models/'));
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
raw_data$wild_type[[paste("results_",10,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));

exp_dirs <- Sys.glob(file.path('../../results/S178A/*/adhesion_props/models/'));
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
raw_data$S178A[[paste("results_",10,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));

print('Done Loading Data')

########################################
#Result filtering
########################################
processed = list();

for (exp_type in names(raw_data)) {
    for (n in names(raw_data[[exp_type]])) {
    	processed$no_filt[[exp_type]][[n]] = filter_results(raw_data[[exp_type]][[n]], 
        	min_R_sq = -Inf, max_p_val = Inf);
	    
	    #now figure out the minimum length of the experiment and check to make
	    #sure that is indeed the min
	    regex_hit = regexpr('[[:digit:]]+$',n);
	    min_length = substr(n, regex_hit[1], regex_hit[1] + attr(regex_hit,"match.length"));

	    stopifnot(min(processed$no_filt[[exp_type]][[n]]$assembly$length) == min_length);
	    stopifnot(min(processed$no_filt[[exp_type]][[n]]$disassembly$length) == min_length);
	
		processed$only_signif[[exp_type]][[n]] = filter_results(raw_data[[exp_type]][[n]], 
			min_R_sq = -Inf, max_p_val = 0.05);
    }
}

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
stop()

################################################################################
#Plotting
################################################################################
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);
phase_lengths = seq(5,10);
ratio_conf_str = list();

#Assembly phase length variation
svg(file.path(out_folder,'controls','assembly_length_variation.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
par(bty='n',mar=c(3,4.3,2,0))

for (length in phase_lengths) {
    results_str = paste("results_", length,sep = '');
    
    data_points = list(processed$only_signif$wild_type[[results_str]]$assembly$slope,
        processed$only_signif$S178A[[results_str]]$assembly$slope);

    ratio_conf_str$assembly[[results_str]] = print_ratio_conf_string(data_points[[1]], data_points[[2]]);

    boxplot_with_points(data_points,
        names=c('Wild-type','S178A'), ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')), 
        main=paste("Minimum Phase Length:",length), inc.points=FALSE)
}
graphics.off()

#Disassembly phase length variation
svg(file.path(out_folder,'controls','disassembly_length_variation.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
par(bty='n',mar=c(3,4.3,2,0))

for (length in phase_lengths) {
    results_str = paste("results_", length,sep = '');

    data_points = list(processed$only_signif$wild_type[[results_str]]$disassembly$slope,
        processed$only_signif$S178A[[results_str]]$disassembly$slope);

    ratio_conf_str$disassembly[[results_str]] = print_ratio_conf_string(data_points[[1]], data_points[[2]]);
    
    boxplot_with_points(data_points,
        names=c('Wild-type','S178A'), ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')), 
        main=paste("Minimum Phase Length:",length), inc.points=FALSE)
}
graphics.off()

#Edge Dist at birth variation
svg(file.path(out_folder,'controls','edge_dist_at_birth.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
par(bty='n',mar=c(3,4.3,2,0))

for (length in phase_lengths) {
    results_str = paste("results_", length,sep = '');
    
    data_points = list(processed$only_signif$wild_type[[results_str]]$assembly$edge_dist,
        processed$only_signif$S178A[[results_str]]$assembly$edge_dist);

    ratio_conf_str$birth_dist[[results_str]] = print_ratio_conf_string(data_points[[1]], data_points[[2]]);

    boxplot_with_points(data_points,
        names=c('Wild-type','S178A'), ylab=expression(paste('Distance from Cell Edge (',mu,'m)',sep='')), 
        main=paste("Minimum Phase Length:",length), inc.points=FALSE)
}
graphics.off()

#Edge dist at death variation
svg(file.path(out_folder,'controls','edge_dist_at_death.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
par(bty='n',mar=c(3,4.3,2,0))

for (length in phase_lengths) {
    results_str = paste("results_", length,sep = '');
    
    data_points = list(processed$only_signif$wild_type[[results_str]]$disassembly$edge_dist,
        processed$only_signif$S178A[[results_str]]$disassembly$edge_dist);

    ratio_conf_str$death_dist[[results_str]] = print_ratio_conf_string(data_points[[1]], data_points[[2]]);

    boxplot_with_points(data_points,
        names=c('Wild-type','S178A'), ylab=expression(paste('Distance from Cell Edge (',mu,'m)',sep='')), 
        main=paste("Minimum Phase Length:",length), inc.points=FALSE)
}
graphics.off()
