rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
debug=TRUE;

################################################################################
#Result loading
################################################################################
raw_data <- list()
single_props <- list()

#Wild-type FA
exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

raw_data$wild_type$intensity = load_results(exp_dirs,file.path('intensity.Rdata'));
raw_data$wild_type$static_props <- load_data_files(exp_dirs, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

raw_data$wild_type$average_cell_intensity = load_data_files(exp_dirs,c('../single_props/Cell_mean_intensity.csv'),inc_exp_names=F)
raw_data$wild_type$average_noad_intensity = load_data_files(exp_dirs,c('../single_props/Cell_not_ad_mean_intensity.csv'),inc_exp_names=F)
raw_data$wild_type$average_ad_intensity = load_data_files(exp_dirs,c('../single_props/Adhesion_mean_intensity.csv'),inc_exp_names=F)

#S178A Results
exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]

raw_data$S178A$intensity = load_results(exp_dirs_S,file.path('intensity.Rdata'));
raw_data$S178A$static_props <- load_data_files(exp_dirs_S, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

raw_data$S178A$average_cell_intensity = load_data_files(exp_dirs_S,c('../single_props/Cell_mean_intensity.csv'),inc_exp_names=F)
raw_data$S178A$average_noad_intensity = load_data_files(exp_dirs_S,c('../single_props/Cell_not_ad_mean_intensity.csv'),inc_exp_names=F)
raw_data$S178A$average_ad_intensity = load_data_files(exp_dirs_S,c('../single_props/Adhesion_mean_intensity.csv'),inc_exp_names=F)

#Wild-type FAK
exp_dirs <- Sys.glob('../../results/FAK/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

raw_data$FAK$intensity = load_results(exp_dirs,file.path('intensity.Rdata'));
raw_data$FAK$static_props <- load_data_files(exp_dirs, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

print('Done Loading Data')

########################################
#Result filtering
########################################

processed = list();
dynamic_props = list();
static_props = list();
for (exp_type in names(raw_data)) {
    for (property in names(raw_data[[exp_type]])) {
        if (property != "intensity") {
            next;
        }
        if (debug) {
            print(paste("Filtering", exp_type, property));
        }

        processed$no_filt[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
            min_R_sq = -Inf, max_p_val = Inf, pos_slope=FALSE);
        
        processed$only_signif[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
            min_R_sq = -Inf, max_p_val = 0.05);
    }
    
    dynamic_props[[exp_type]] = gather_general_dynamic_props(raw_data[[exp_type]]$intensity)
    static_props[[exp_type]] = gather_static_props(raw_data[[exp_type]]$static_props) 
}

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

find_mean_from_list <- function(data) {
    this_mean_set = lapply(data,as.numeric);
    this_mean_set = unlist(lapply(this_mean_set,mean));
    return(this_mean_set);
}

build_column_matrix <- function(data) {
    max_size <- max(unlist(lapply(data,length)));

    data_mat = matrix(NA,nrow=max_size,ncol=length(data))

    for (i in 1:length(data)) {
        data_mat[1:length(data[[i]]),i] = t(data[[i]])
    }

    return(data_mat);
}

stop()

################################################################################
#Plotting
################################################################################

wt_exp_average_int = find_mean_from_list(raw_data$wild_type$average_cell_intensity);
wt_exp_average_noad_int = find_mean_from_list(raw_data$wild_type$average_noad_intensity);
wt_exp_average_ad_int = find_mean_from_list(raw_data$wild_type$average_ad_intensity);

S178A_exp_average_int = find_mean_from_list(raw_data$S178A$average_cell_intensity);
S178A_exp_average_noad_int = find_mean_from_list(raw_data$S178A$average_noad_intensity);
S178A_exp_average_ad_int = find_mean_from_list(raw_data$S178A$average_ad_intensity);

dir.create(dirname(file.path(out_folder,'per_cell_comparisons','fluor_intensity.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'per_cell_comparisons','fluor_intensity.svg'), height=7/2, width=7*(3/2));
layout(rbind(c(1,2,3)))
par(bty='n',mar=c(2.2,4,0.6,0))

boxplot_with_points(list(wt_exp_average_int, S178A_exp_average_int),names=c('Wild-type','S178A'),
    ylab='Average Cellular Fluorescence')
boxplot_with_points(list(wt_exp_average_noad_int, S178A_exp_average_noad_int),names=c('Wild-type','S178A'),
    ylab='Average Non-adhesion Cell Fluorescence')
boxplot_with_points(list(wt_exp_average_ad_int, S178A_exp_average_ad_int),names=c('Wild-type','S178A'),
    ylab='Average Adhesion Fluorescence')

graphics.off()

wt_full_cell_mat  = build_column_matrix(raw_data$wild_type$average_cell_intensity)
S178A_full_cell_mat  = build_column_matrix(raw_data$S178A$average_cell_intensity)

svg(file.path(out_folder,'per_cell_comparisons','cell_intensity_seq.svg'), height=7/2, width=7);
layout(rbind(c(1,2)))
par(bty='n',mar=c(4,4,0,0))

matplot(wt_full_cell_mat,type='l', xlab='Time (minutes)', ylab='Average Cellular Fluorescence')

par(bty='n',mar=c(4.2,4,0,0))
matplot(S178A_full_cell_mat,type='l', xlab='Time (minutes)', ylab='Average Cellular Fluorescence')

graphics.off()
