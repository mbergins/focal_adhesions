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

exp_dirs <- Sys.glob('../../results/FA_RhoG/*')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

for (exp_type_dir in exp_dirs) {
    if (basename(exp_type_dir) == "KD-Problems") {
        next;
    }
    exp_name = basename(exp_type_dir)
    model_folder_set = Sys.glob(file.path(exp_type_dir,'*/*/models'))
    raw_data[[exp_name]]$intensity = load_results(model_folder_set,file.path('intensity.Rdata'), debug=TRUE)
    raw_data[[exp_name]]$static_props <- load_data_files(model_folder_set, 
        file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);
}
print('Done Loading Data')
########################################
#Result filtering
########################################
processed = list();
dynamic_props = list();
static_props = list();

for (exp_name in names(raw_data)) {
    # for (property in names(raw_data[[exp_type]])) {
    #     if (property == "static_props") {
    #         next;
    #     }
    #     if (debug) {
    #         print(paste("Filtering", exp_name));
    #     }
    # }
    processed$no_filt[[exp_name]] = filter_results(raw_data[[exp_name]]$intensity, 
        min_R_sq = -Inf, max_p_val = Inf);
    
    processed$only_signif[[exp_name]] = filter_results(raw_data[[exp_name]]$intensity, 
        min_R_sq = -Inf, max_p_val = 0.05);

    dynamic_props[[exp_name]] = gather_general_dynamic_props(raw_data[[exp_name]]$intensity)
    static_props[[exp_name]] = gather_static_props(raw_data[[exp_name]]$static_props) 
}

#rm(raw_data)
gc()

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);
stop();
################################################################################
#Plotting
################################################################################
dir.create(dirname(file.path(out_folder, 'RhoG', 'kinetics_comparison.svg')), 
    recursive=TRUE, showWarnings=FALSE);

stage_data <- gather_stage_lengths(processed$only_signif$control, 
                                   processed$only_signif$KD)

#Static properties
svg(file.path(out_folder, 'RhoG', 'statics_comparison.svg'));
par(mar=c(3,4.5,1,0), bty='n');
boxplot(static_props$control$Area, static_props$KD$Area, names=c('Control', 'KD'), ylab = 'FA Area (\u03BCm\u00B2)', range=0)

t_test_results = t.test(static_props$control$Area, static_props$KD$Area);

if (t_test_results$p.value == 0) {
    t_test_results$p.value = 1e-05;
}

plot_dims = par("usr");
x_pos = (plot_dims[2] - plot_dims[1])*0.5 + plot_dims[1]
y_pos = (plot_dims[4] - plot_dims[3])*0.9 + plot_dims[3]

text(x_pos,y_pos, paste('p<',t_test_results$p.value, sep=''), col="blue");

graphics.off()

#made this file for Chris W. so he can do his own plotting (2/8/2010)
# area_temp = data.frame(control = c(static_props$control$Area, rep(NA, 79147-73992)), KD = static_props$KD$Area);
# write.csv(area_temp, file="RhoG_KD_Areas.csv", row.names=F);

#Static Property Barplot
control_area_props = t.test(static_props$control$Area);
KD_area_props = t.test(static_props$KD$Area);

max_conf_int = max(c(control_area_props$conf.int, KD_area_props$conf.int))*1.01;

svg(file.path(out_folder, 'RhoG', 'statics_barplot.svg'))
par(mar=c(2,4,0,0))
barplot(c(control_area_props$estimate, KD_area_props$estimate), 
    ylim=c(0,max_conf_int), names=c('Control', 'KD'), ylab="Mean Adhesion Area (\u03BCm\u00B2)");
errbar(c(0.7,1.9), c(control_area_props$estimate, KD_area_props$estimate), 
    c(control_area_props$conf.int[1], KD_area_props$conf.int[1]),
    c(control_area_props$conf.int[2], KD_area_props$conf.int[2]), add=TRUE, cap=0.25, cex=1e-10
)
graphics.off();


#Comparing the assembly/disassembly rates
dir.create(dirname(file.path(out_folder, 'RhoG', 'kinetics_comparison.svg')), 
    recursive=TRUE, showWarnings=FALSE);

svg(file.path(out_folder, 'RhoG', 'kinetics_comparison.svg'), height=7*(1/2));
layout(rbind(c(1,2)))
par(mar=c(4,4.5,1,0), bty='n');
boxplot_with_points(list(processed$only_signif$KD$assembly$slope, 
                         processed$only_signif$control$assembly$slope), 
                    names=c("KD", "Control"), inc.points=FALSE,
                    ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
boxplot_with_points(list(processed$only_signif$KD$disassembly$slope, 
                         processed$only_signif$control$disassembly$slope), 
                    names=c("KD", "Control"), inc.points=FALSE,
                    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
graphics.off()

# for (i in unique(processed$only_signif$KD$assembly$exp_num)) {
#     for (j in unique(processed$only_signif$KD$assembly$exp_num)) {
#         if (j == i) {
#             next;
#         }
#         if (i > j) {
#             next;
#         }
#         filtered_data = subset(processed$only_signif$KD$assembly, exp_num != i | exp_num != j);
#         filtered_dis = subset(processed$only_signif$KD$disassembly, exp_num != i | exp_num != j);
#         
#         print(paste(i,j))
# 
#         p_vals = determine_median_p_value(filtered_data$slope, processed$only_signif$control$assembly$slope);
#         print(p_vals$ratio_conf)
#         p_vals = determine_median_p_value(filtered_dis$slope, processed$only_signif$control$disassembly$slope);
#         print(p_vals$ratio_conf)
#     }
# }

# svg(file.path(out_folder, 'RhoG', 'single_exp.svg'));
# layout(rbind(c(1,2),c(3,4)))
# par(mar=c(4,4.5,1,0), bty='n');
# kd_exp_assem = boxplot(slope ~ exp_num, data=processed$only_signif$KD$assembly, 
#     notch=T, varwidth=T,ylim=c(0,0.12), xlab='KD Experiments',
#     ylab=expression(paste('Assembly Rate (',min^-1,')')), pch=19,cex=0.25)
# boxplot(slope ~ exp_num, data=processed$only_signif$control$assembly, 
#     notch=T, varwidth=T, ylim=c(0,0.12), xlab='Control Experiments',
#     ylab=expression(paste('Assembly Rate (',min^-1,')')), pch=19,cex=0.25)
# 
# kd_exp_disassem = boxplot(slope ~ exp_num, data=processed$only_signif$KD$disassembly, 
#     notch=T, varwidth=T, ylim=c(0,0.12), xlab='KD Experiments',
#     ylab=expression(paste('Disassembly Rate (',min^-1,')')), pch=19,cex=0.25)
# boxplot(slope ~ exp_num, data=processed$only_signif$control$disassembly, 
#     notch=T, varwidth=T, ylim=c(0,0.12),xlab='Control Experiments',
#     ylab=expression(paste('Disassembly Rate (',min^-1,')')), pch=19,cex=0.25)
# 
# graphics.off()
# 
