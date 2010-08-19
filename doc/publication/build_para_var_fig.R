rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
debug = TRUE;

################################################################################
#Result loading
################################################################################
raw_data <- list()

exp_dirs <- Sys.glob('../../results/ID_thresh/*')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

for (exp_type_dir in exp_dirs) {
    exp_name = basename(exp_type_dir)
    thresh_levels = Sys.glob(file.path(exp_type_dir,'*'))
    for (dir_set in thresh_levels) {
        thresh_level = basename(dir_set)
        model_folder_set = Sys.glob(file.path(dir_set,'*/*/models'))
        raw_data[[exp_name]][[thresh_level]] = load_results(model_folder_set,file.path('intensity.Rdata'), debug=TRUE)
    }
}

reg_thresh_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/models/')
reg_thresh_dirs <- reg_thresh_dirs[file_test('-d',reg_thresh_dirs)]

raw_data[["FA"]][["0_10"]] = load_results(reg_thresh_dirs,file.path('intensity.Rdata'), debug=TRUE)

exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]

raw_data[["S178A"]][["0_10"]] = load_results(exp_dirs_S,file.path('intensity.Rdata'), debug=TRUE)


exp_dirs_reduced <- Sys.glob('../../results/focal_adhesions_reduced/*/adhesion_props/models/')
exp_dirs_reduced <- exp_dirs_reduced[file_test('-d',exp_dirs_reduced)]
raw_data[["FA_reduced"]] = load_results(exp_dirs_reduced,file.path('intensity.Rdata'))

exp_dirs_S <- Sys.glob('../../results/S178A_reduced/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]
raw_data[["S178A_reduced"]] = load_results(exp_dirs_S,file.path('intensity.Rdata'))

exp_dirs <- Sys.glob('../../results/focal_adhesions_no_split/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data[["FA_no_split"]] = load_results(exp_dirs,file.path('intensity.Rdata'))

exp_dirs <- Sys.glob('../../results/S178A_no_split/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data[["S178A_no_split"]] = load_results(exp_dirs,file.path('intensity.Rdata'))

print('Done Loading Data')

########################################
#Result filtering
########################################
processed = list();

for (exp_name in names(raw_data)) {
    if (any(exp_name == c("FA_reduced", "S178A_reduced", "FA_no_split", "S178A_no_split"))) {
        processed$no_filt[[exp_name]] = filter_results(raw_data[[exp_name]], 
            min_R_sq = -Inf, max_p_val = Inf);
        
        processed$only_signif[[exp_name]] = filter_results(raw_data[[exp_name]], 
            min_R_sq = -Inf, max_p_val = 0.05);
    }
    for (thresh_level in names(raw_data[[exp_name]])) {
        if (debug) {
            print(paste("Filtering", exp_name, thresh_level));
        }
        processed$no_filt[[exp_name]][[thresh_level]] = filter_results(raw_data[[exp_name]][[thresh_level]], 
            min_R_sq = -Inf, max_p_val = Inf);
        
        processed$only_signif[[exp_name]][[thresh_level]] = filter_results(raw_data[[exp_name]][[thresh_level]], 
            min_R_sq = -Inf, max_p_val = 0.05);
    }
}

#rm(raw_data)
gc()

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);
stop()

################################################################################
#Plotting
################################################################################

########################################
# FA threshold variation
########################################
dir.create(file.path(out_folder, 'parameter_variation'), recursive=TRUE, showWarnings=FALSE);

svg(file.path(out_folder, 'parameter_variation', 'comparison_assembly.svg'), height=7*(3/2));
layout(rbind(c(1,1,2,2),c(3,3,4,4), c(0,5,5,0)))
par(mar=c(4,4.5,1,0), bty='n');

thresh_names = c("0_06","0_07","0_08","0_09","0_10");
thresh_decimals = c("0.06","0.07","0.08","0.09","0.10");
for (i in 1:length(thresh_names)) {
    boxplot_with_points(list(processed$only_signif$FA[[thresh_names[i]]]$assembly$slope, 
                             processed$only_signif$S178A[[thresh_names[i]]]$assembly$slope), 
                        names=c('wild-type','S178A'), inc.points=F, main=paste('Threshold',thresh_decimals[i]),
                        ylab=expression(paste('Assembly rate (',min^-1,')',sep='')), yaxt='n', ylim=c(0,0.16)); 
    axis(2,at=seq(0,0.16,by=0.02));
}
graphics.off()

print("Done with threshold variation vs assembly rate")

svg(file.path(out_folder, 'parameter_variation', 'comparison_disassembly.svg'), height=7*(3/2))
layout(rbind(c(1,1,2,2),c(3,3,4,4), c(0,5,5,0)))
par(mar=c(4,4.5,1,0), bty='n');

for (i in 1:length(thresh_names)) {
    boxplot_with_points(list(processed$only_signif$FA[[thresh_names[i]]]$disassembly$slope, 
                             processed$only_signif$S178A[[thresh_names[i]]]$disassembly$slope), 
                        names=c('wild-type','S178A'), inc.points=F, main=paste('Threshold',thresh_decimals[i]),
                        ylab=expression(paste('Disassembly rate (',min^-1,')',sep='')), yaxt='n', ylim=c(0,0.16)); 
    axis(2,at=seq(0,0.14,by=0.02));
}

graphics.off()
print("done with threshold variation vs disassembly rate")

# svg(file.path(out_folder, 'thresh_variation', 'wild_type_assembly.svg'), height=7, width=7*1/2)
# layout(rbind(1,2))
# boxplot_with_points(list(processed$only_signif$FA[["0_06"]]$assembly$slope, 
#                          processed$only_signif$FA[["0_07"]]$assembly$slope,
#                          processed$only_signif$FA[["0_08"]]$assembly$slope,
#                          processed$only_signif$FA[["0_09"]]$assembly$slope,
#                          processed$only_signif$FA[["0_10"]]$assembly$slope), 
#                     names=seq(0.06,0.1,length=5), with.median.props=FALSE, notch=T)
# boxplot_with_points(list(processed$only_signif$FA[["0_06"]]$assembly$slope, 
#                          processed$only_signif$FA[["0_07"]]$assembly$slope,
#                          processed$only_signif$FA[["0_08"]]$assembly$slope,
#                          processed$only_signif$FA[["0_09"]]$assembly$slope,
#                          processed$only_signif$FA[["0_10"]]$assembly$slope), 
#                     names=seq(0.06,0.1,length=5), with.median.props=FALSE, notch=T)
# 
# boxplot_with_points(list(processed$only_signif$S178A[["0_06"]]$assembly$slope, 
#                          processed$only_signif$S178A[["0_07"]]$assembly$slope,
#                          processed$only_signif$S178A[["0_08"]]$assembly$slope,
#                          processed$only_signif$S178A[["0_09"]]$assembly$slope,
#                          processed$only_signif$S178A[["0_10"]]$assembly$slope), 
#                     names=seq(0.06,0.1,length=5), with.median.props=FALSE, notch=T)


########################################
# Image resampling results
########################################

#Reduced Wild-type
svg(file.path(out_folder,'parameter_variation','resampled_wt_slopes.svg'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(3,4.2,0,0));

boxplot_with_points(list(processed$only_signif$FA$"0_10"$assembly$slope,
                         processed$only_signif$FA_reduced$assembly$slope/2), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')), 
                    inc.points=F)

boxplot_with_points(list(processed$only_signif$FA$"0_10"$disassembly$slope,
                         processed$only_signif$FA_reduced$disassembly$slope/2), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')), 
                    inc.points=F)

boxplot_with_points(list(processed$only_signif$FA$"0_10"$assembly$R_sq,
                         processed$only_signif$FA_reduced$assembly$R_sq), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Assembly ',R^2,sep='')), 
                    inc.points=F, median.props.pos=c(0.5,0.2))

boxplot_with_points(list(processed$only_signif$FA$"0_10"$disassembly$R_sq,
                         processed$only_signif$FA_reduced$disassembly$R_sq), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Disassembly ',R^2,sep='')), 
                    inc.points=F, median.props.pos=c(0.5,0.2))
graphics.off()
print("done with image resampling in wild-type")

#Reduced S178A
svg(file.path(out_folder,'parameter_variation','resampled_S178A_slopes.svg'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(3,4.2,0,0));

boxplot_with_points(list(processed$only_signif$S178A$"0_10"$assembly$slope,
                         processed$only_signif$S178A_reduced$assembly$slope/2), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
                    inc.points=F)

boxplot_with_points(list(processed$only_signif$S178A$"0_10"$disassembly$slope,
                         processed$only_signif$S178A_reduced$disassembly$slope/2), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
                    inc.points=F)

boxplot_with_points(list(processed$only_signif$S178A$"0_10"$assembly$R_sq,
                         processed$only_signif$S178A_reduced$assembly$R_sq), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Assembly Adjusted ',R^2,sep='')), 
                    median.props.pos=c(0.5,0.2), inc.points=F)

boxplot_with_points(list(processed$only_signif$S178A$"0_10"$disassembly$R_sq,
                         processed$only_signif$S178A_reduced$disassembly$R_sq), 
                    names=c('All', 'Sampled'), ylab=expression(paste('Disassembly Adjusted ',R^2,sep='')), 
                    median.props.pos=c(0.5,0.2), inc.points=F)
graphics.off()
print("done with image resampling in S178A")
