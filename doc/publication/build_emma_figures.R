rm(list = ls())
start_time = proc.time()
library(lattice)
library(geneplotter)
library(Hmisc)
source('FA_analysis_lib.R')
source('bilinear_modeling.R')
source('FA_alignment_search.R')

################################################################################
#Result loading
################################################################################
raw_data <- list()
raw_data_base_dir = '../../results/emma/processed_2stdev_crop/';

########################################
# Fibronectin Percentages - NS
########################################

#Fibro 1ug - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_trial_2/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Fibro_1ug$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#Fibro 10ug - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_trial_2/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Fibro_10ug$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#Fibro 100ug - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_trial_2/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Fibro_100ug$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

########################################
# Fibronectin Percentages - 2xKD
########################################

#Fibro 1ug - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_2xKD/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Fibro_1ug_2xKD$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#Fibro 10ug - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_2xKD/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Fibro_10ug_2xKD$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#Fibro 100ug - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_2xKD/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Fibro_100ug_2xKD$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

########################################
# Rat2 
########################################

#Rat2 CK-689
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Rat2/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Rat2$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#Rat2 CK-666
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'Rat2_ck*/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$Rat2_ck$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

print('Done Loading Data')

###########################################################
#Result filtering
###########################################################
source('FA_analysis_lib.R')
dynamic_props = list();
adhesion_counts = list();
for (exp_type in names(raw_data)) {
    dynamic_props[[exp_type]] = gather_general_dynamic_props(raw_data[[exp_type]]$intensity, min.longevity=1,
        debug=F)
    adhesion_counts[[exp_type]] = count_adhesions_per_image(raw_data[[exp_type]]$intensity,time.spacing=2.5)
}


print('Done Filtering Data')
out_folder = '../../doc/publication/figures/emma'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);
 
system("notify-send \"done reading in R data\"");
stop()

################################################################################
#Plotting
################################################################################

source('FA_analysis_lib.R')

png(file.path(out_folder,'FA_summary','area.png'));
par(bty='n',mar=c(2.8,2.6,.5,0), mgp=c(1.6,0.5,0),xpd=T)
plot_dual_hist(dynamic_props$Fibro_1ug$mean_area,dynamic_props$Fibro_100ug$mean_area,break.int=0.25)
axis(1)
axis(2)
title(xlab='Mean Adhesion Area',ylab='Density')
legend('topright',c('Fibro 1ug','Fibro 100ug'),fill=c(rgb(0,0,1,0.6),rgb(1,1,0,0.6)),bty='n')
graphics.off()

png(file.path(out_folder,'FA_summary','axial_ratio.png'));
par(bty='n',mar=c(2.8,2.6,.5,0), mgp=c(1.6,0.5,0),xpd=T)
plot_dual_hist(dynamic_props$Fibro_1ug$mean_axial_ratio,dynamic_props$Fibro_100ug$mean_axial_ratio)
axis(1)
axis(2)
title(xlab='Mean Axial Ratio',ylab='Density')
legend('topright',c('Fibro 1ug','Fibro 100ug'),fill=c(rgb(0,0,1,0.6),rgb(1,1,0,0.6)),bty='n')
graphics.off()

png(file.path(out_folder,'FA_summary','longevity.png'));
par(bty='n',mar=c(2.8,2.6,.5,0), mgp=c(1.6,0.5,0),xpd=T)
plot_dual_hist(dynamic_props$Fibro_1ug$longevity,dynamic_props$Fibro_100ug$longevity,break.int=5)
axis(1)
axis(2)
title(xlab='Longevity (min)',ylab='Density')
legend('topright',c('Fibro 1ug','Fibro 100ug'),fill=c(rgb(0,0,1,0.6),rgb(1,1,0,0.6)),bty='n')
graphics.off()

###########################################################
# Global FA Properties
###########################################################

dynamic_props_names = names(dynamic_props$Fibro_1ug);
prop_summaries = list()

source('FA_analysis_lib.R')
for (prop in dynamic_props_names) {
    default.table.text.format = '%.1f';
    if (prop == "mean_area") {
        default.table.text.format = '%.3f';
    }
    # default.table.text.format = '%.1f\u00B1%.2f';
    # if (prop == "mean_area") {
    #     default.table.text.format = '%.3f\u00B1%.5f';
    # }
    # if (prop == "average_speed") {
    #     default.table.text.format = '%.3f\u00B1%.4f';
    # }
    # if (prop == "mean_axial_ratio") {
    #     default.table.text.format = '%.1f\u00B1%.4f';
    # }
    # if (prop == "mean_major_axis") {
    #     default.table.text.format = '%.1f\u00B1%.4f';
    # }
    # if (prop == "mean_minor_axis") {
    #     default.table.text.format = '%.1f\u00B1%.4f';
    # }
    
    data_set = list(dynamic_props$Fibro_1ug[[prop]],
        dynamic_props$Fibro_10ug[[prop]],dynamic_props$Fibro_100ug[[prop]],
        dynamic_props$Fibro_1ug_2xKD[[prop]],dynamic_props$Fibro_10ug_2xKD[[prop]],
        dynamic_props$Fibro_100ug_2xKD[[prop]],dynamic_props$Rat2[[prop]],
        dynamic_props$Rat2_ck[[prop]]);

    prop_summaries[[prop]] = gather_dynamics_summary_median(data_set,
        table.text.format = default.table.text.format);
}
prop_summaries$exp_num = lapply(dynamic_props, function (x) length(unique(x$exp_num)))

ad_count_summary = gather_dynamics_summary(list(adhesion_counts$Fibro_1ug$counts,
        adhesion_counts$Fibro_10ug$counts,adhesion_counts$Fibro_100ug$counts,
        adhesion_counts$Fibro_1ug_2xKD$counts,adhesion_counts$Fibro_10ug_2xKD$counts,
        adhesion_counts$Fibro_100ug_2xKD$counts, adhesion_counts$Rat2$counts,
        adhesion_counts$Rat2_ck$counts),table.text.format='%.1f\u00B1%.1f');

props_of_interest = c('mean_area','mean_axial_ratio','mean_major_axis','mean_minor_axis',
    'longevity')
# for (prop in props_of_interest) {
#     data = prop_summaries[[prop]];
#     
#     svg(file.path(out_folder,'overall_FA_props',paste(prop,'.svg',sep='')),width=4,height=5)
#     par(bty='n',mar=c(2.8,2.6,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
# 
#     bar_data = barplot(unlist(data$mean), ylim=c(0,max(unlist(data$t_conf))),
#         ylab=prop,names=c('1','10','100','1','10','100','Rat2','Rat2 CK'))
# 
#     for (i in 1:length(data$mean)) {
#         errbar(bar_data[i],data$mean[[i]],data$t_conf[[i]][2],
#             data$t_conf[[i]][1],add=T,cex=1E-10)
#     }
# 
#     graphics.off()
# }

table_data = list()
for (prop in props_of_interest) {
    table_data[[prop]] = unlist(prop_summaries[[prop]]$table_text)
}
table_data$ad_per_cell = unlist(ad_count_summary$table_text)
table_data$counts = as.character(unlist(prop_summaries$longevity$counts))
table_data$exp_num = as.character(prop_summaries$exp_num)

table_data_frame = data.frame(table_data,
    row.names=c('NS 1\u03BCg/mL','NS 10\u03BCg/mL','NS 100\u03BCg/mL',
        '2xKD 1\u03BCg/mL','2xKD 10\u03BCg/mL','2xKD 100\u03BCg/mL',
        'Rat2 CK-689','Rat2 CK-666'))
table_data_frame = t(table_data_frame);
row.names(table_data_frame) <- c('Mean Area (\u03BCm)','Mean Axial Ratio','Mean Major Axis (\u03BCm)',
    'Mean Minor Axis (\u03BCm)','Mean Longevity (min)', 'Adhesions/Cell/10 min','# Adhesions','# Cells')

write.table(table_data_frame,file='../../doc/publication/figures/emma/FA_summary/FA_summary_median.csv',sep=',')
