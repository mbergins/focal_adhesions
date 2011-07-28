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
raw_data_base_dir = '../../results/emma/processed_2stdev/';

########################################
# Various adhesion tags, all on 100ug fibronectin
########################################

#KD - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_KD/Pax_*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$KD_Pax$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#IA32 - Pax
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/Pax_*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$WT_Pax$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#KD - FAK 
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_KD/FAK_*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$KD_FAK$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#IA32 - FAK 
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/FAK_*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$WT_FAK$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#KD - Vin
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_KD/Vin_*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$KD_Vin$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

#IA32 - Vin
exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/Vin_*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$WT_Vin$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

########################################
# Fibronectin Percentages
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

print('Done Loading Data')

###########################################################
#Result filtering
###########################################################
total_models = 0;
for (exp_type in names(raw_data)) {
    total_models = total_models + find_number_of_models(raw_data[[exp_type]]$intensity);
}

processed = list();
dynamic_props = list();
dynamic_props_five = list();
full_cell_props = list();
static_props = list();
for (exp_type in names(raw_data)) {
    # for (property in names(raw_data[[exp_type]])) {
    #     if (property != "intensity") {
    #         next;
    #     }
    #     print(paste("Filtering", exp_type, property));
    #     
    #     processed$no_filt[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
    #         min.r.sq = -Inf, max.p.val = Inf, pos.slope=FALSE);
    #     
    #     processed$only_signif[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
    #         min.r.sq = -Inf, max.p.val = 0.05, model_count = total_models);
    # }
    
    # full_cell_props[[exp_type]] = gather_global_exp_summary(raw_data[[exp_type]]$intensity)
    dynamic_props[[exp_type]] = gather_general_dynamic_props(raw_data[[exp_type]]$intensity, min.longevity=NA,
        debug=F)
    dynamic_props_five[[exp_type]] = gather_general_dynamic_props(raw_data[[exp_type]]$intensity, 
        min.longevity=5, debug=F)
}

print('Done Filtering Data')
out_folder = '../../doc/publication/figures/emma'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

KD_Pax = processed$only_signif$KD_Pax$intensity;
IA32_Pax = processed$only_signif$WT_Pax$intensity;

KD_Vin = processed$only_signif$KD_Vin$intensity;
IA32_Vin = processed$only_signif$WT_Vin$intensity;

KD_FAK = processed$only_signif$KD_FAK$intensity;
IA32_FAK = processed$only_signif$WT_FAK$intensity;

Fibro_1ug = processed$only_signif$Fibro_1ug$intensity;
Fibro_10ug = processed$only_signif$Fibro_10ug$intensity;
Fibro_100ug = processed$only_signif$Fibro_100ug$intensity;

# rate_angle_data = list()
# rate_angle_data$Fibro_100ug = get_angle_data_set(raw_data$Fibro_100ug$intensity)
# rate_angle_data$Fibro_10ug = get_angle_data_set(raw_data$Fibro_10ug$intensity)
# rate_angle_data$Fibro_1ug = get_angle_data_set(raw_data$Fibro_1ug$intensity)
# 
# split_data = list()
# split_data$Fibro_100ug = split_angle_data(rate_angle_data$Fibro_100ug)
# split_data$Fibro_10ug = split_angle_data(rate_angle_data$Fibro_10ug)
# split_data$Fibro_1ug = split_angle_data(rate_angle_data$Fibro_1ug)
 
stop()
 
###########################################################
# P-value calculations
###########################################################

Pax_p = list();
Pax_p$assembly = determine_median_p_value(IA32_Pax$assembly$slope,KD_Pax$assembly$slope)
Pax_p$disassembly = determine_median_p_value(IA32_Pax$disassembly$slope,KD_Pax$disassembly$slope)

FAK_p = list();
FAK_p$assembly = determine_median_p_value(IA32_FAK$assembly$slope,KD_FAK$assembly$slope)
FAK_p$disassembly = determine_median_p_value(IA32_FAK$disassembly$slope,KD_FAK$disassembly$slope)

Vin_p = list();
Vin_p$assembly = determine_median_p_value(IA32_Vin$assembly$slope,KD_Vin$assembly$slope)
Vin_p$disassembly = determine_median_p_value(IA32_Vin$disassembly$slope,KD_Vin$disassembly$slope)

fibro_1ug_10ug_p = list()
fibro_1ug_10ug_p$assembly = determine_median_p_value(Fibro_1ug$assembly$slope,Fibro_10ug$assembly$slope)
fibro_1ug_10ug_p$disassembly = determine_median_p_value(Fibro_1ug$disassembly$slope,Fibro_10ug$disassembly$slope)

fibro_1ug_100ug_p = list()
fibro_1ug_100ug_p$assembly = determine_median_p_value(Fibro_1ug$assembly$slope,Fibro_100ug$assembly$slope)
fibro_1ug_100ug_p$disassembly = determine_median_p_value(Fibro_1ug$disassembly$slope,Fibro_100ug$disassembly$slope)

fibro_10ug_100ug_p = list()
fibro_10ug_100ug_p$assembly = determine_median_p_value(Fibro_10ug$assembly$slope,Fibro_100ug$assembly$slope)
fibro_10ug_100ug_p$disassembly = determine_median_p_value(Fibro_10ug$disassembly$slope,Fibro_100ug$disassembly$slope)

prep_time = proc.time()
system("notify-send \"done reading in R data\"");
stop()

################################################################################
#Plotting
################################################################################

########################################
# Average Adhesion Area - Fibronectin
########################################
ad_areas = list(dynamic_props$Fibro_1ug$mean_area, dynamic_props$Fibro_10ug$mean_area,
    dynamic_props$Fibro_100ug$mean_area)

fibro_areas = gather_dynamics_summary(ad_areas)

svg(file.path(out_folder,'overall_FA_props','fibro_concens_areas.svg'),width=4,height=5)
par(bty='n',mar=c(2.8,2.6,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
# bar_data = barplot(unlist(ad_areas_mean), names=rep(c('NS', '2xKD'),3))
bar_data = barplot(unlist(fibro_areas$mean), ylim=c(0,max(unlist(fibro_areas$t_conf))),
    ylab='Mean Adhesion Area (pixels)', names=c(1,10,100))

for (i in 1:length(ad_conf_int)) {
    errbar(bar_data[i],fibro_areas$mean[[i]],fibro_areas$t_conf[[i]][2],
        fibro_areas$t_conf[[i]][1],add=T,cex=1E-10)
}

add_labels_with_sub(c(NA,NA,NA),names=c(1,10,100),bar_data[1:3],
    subtitle='Fibronectin Concentration (\u03BCg/mL)',with.axis=F)
graphics.off()

########################################
# Average Axial Ratio - Fibronectin
########################################
fibro_axial = gather_dynamics_summary(list(dynamic_props$Fibro_1ug$mean_axial, 
    dynamic_props$Fibro_10ug$mean_axial,dynamic_props$Fibro_100ug$mean_axial))

svg(file.path(out_folder,'overall_FA_props','fibro_concens_ax_ratio.svg'),width=4,height=5)
par(bty='n',mar=c(2.8,2.6,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
# bar_data = barplot(unlist(ad_areas_mean), names=rep(c('NS', '2xKD'),3))
bar_data = barplot(unlist(fibro_axial$mean), ylim=c(0,max(unlist(fibro_axial$t_conf))),
    ylab='Mean Adhesion Area (pixels)', names=c(1,10,100))

for (i in 1:length(ad_conf_int)) {
    errbar(bar_data[i],fibro_axial$mean[[i]],fibro_axial$t_conf[[i]][2],
        fibro_axial$t_conf[[i]][1],add=T,cex=1E-10)
}

add_labels_with_sub(c(NA,NA,NA),names=c(1,10,100),bar_data[1:3],
    subtitle='Fibronectin Concentration (\u03BCg/mL)',with.axis=F)
graphics.off()


########################################
# Average Adhesion Area - Different Tags
########################################
ad_areas = list(dynamic_props$WT_Pax$mean_area, dynamic_props$KD_Pax$mean_area,
    dynamic_props$WT_Vin$mean_area,dynamic_props$KD_Vin$mean_area,
    dynamic_props$WT_FAK$mean_area,dynamic_props$KD_FAK$mean_area)

tag_areas = gather_dynamics_summary(ad_areas)

svg(file.path(out_folder,'overall_FA_props','tagged_ads_100ug_areas.svg'),width=5,height=5)
par(bty='n',mar=c(2.6,2.6,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
# bar_data = barplot(unlist(ad_areas_mean), names=rep(c('NS', '2xKD'),3))
bar_data = barplot(unlist(tag_areas$mean), ylim=c(0,max(unlist(tag_areas$t_conf))),
    ylab='Mean Adhesion Area (pixels)', names=rep(c('NS', '2xKD'),3))

for (i in 1:length(ad_conf_int)) {
    errbar(bar_data[i],tag_areas$mean[[i]],tag_areas$t_conf[[i]][2],tag_areas$t_conf[[i]][1],add=T,cex=1E-10)
}

add_labels_with_sub(c(NA,NA),names=c('NS','2xKD'),bar_data[1:2],subtitle='Pax',with.axis=F)
add_labels_with_sub(c(NA,NA),names=c('NS','2xKD'),bar_data[3:4],subtitle='Vin',with.axis=F)
add_labels_with_sub(c(NA,NA),names=c('NS','2xKD'),bar_data[5:6],subtitle='FAK',with.axis=F)
graphics.off()

########################################
# Average Axial Ratio - Different Tags
########################################
ad_axials = gather_dynamics_summary(list(dynamic_props$WT_Pax$mean_axial_ratio, 
    dynamic_props$KD_Pax$mean_axial_ratio, dynamic_props$WT_Vin$mean_axial_ratio,
    dynamic_props$KD_Vin$mean_axial_ratio, dynamic_props$WT_FAK$mean_axial_ratio,
    dynamic_props$KD_FAK$mean_axial_ratio))

svg(file.path(out_folder,'overall_FA_props','tagged_ads_100ug_ax_ratio.svg'),width=5,height=5)
par(bty='n',mar=c(2.6,2.6,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
bar_data = barplot(unlist(ad_axials$mean), ylim=c(0,max(unlist(ad_axials$t_conf))),
    ylab='Mean Adhesion Axial Ratio', names=rep(c('NS', '2xKD'),3))

for (i in 1:length(ad_conf_int)) {
    errbar(bar_data[i],ad_axials$mean[[i]],ad_axials$t_conf[[i]][2],ad_axials$t_conf[[i]][1],add=T,cex=1E-10)
}

add_labels_with_sub(c(NA,NA),names=c('NS','2xKD'),bar_data[1:2],subtitle='Pax',with.axis=F)
add_labels_with_sub(c(NA,NA),names=c('NS','2xKD'),bar_data[3:4],subtitle='Vin',with.axis=F)
add_labels_with_sub(c(NA,NA),names=c('NS','2xKD'),bar_data[5:6],subtitle='FAK',with.axis=F)
graphics.off()


########################################
# All-by-All Plots
#######################################

# library(car)
# data_of_interest = c('slope','length','largest_area','ad_sig','mean_axial_ratio','adj.r.squared')
# 
# data_and_file = list('IA32_Pax_assembly.pdf'=IA32_Pax$assembly, 'IA32_Pax_disassembly.pdf'=IA32_Pax$dis,
#     'IA32_Vin_assembly.pdf'=IA32_Vin$assembly, 'IA32_Vin_disassembly.pdf'=IA32_Vin$dis,
#     'KD_Pax_assembly.pdf'=KD_Pax$assembly, 'KD_Pax_disassembly.pdf'=KD_Pax$dis,
#     'KD_Vin_assembly.pdf'=KD_Vin$assembly, 'KD_Vin_disassembly.pdf'=KD_Vin$dis
# )
# 
# for (file_name in names(data_and_file)) {
#     pdf(file.path(out_folder,file_name))
#     lm_sum = summary(lm(slope ~ length + largest_area + ad_sig + mean_axial_ratio + adj.r.squared, 
#         data=data_and_file[[file_name]]))
#     scatterplot.matrix(subset(data_and_file[[file_name]],select = data_of_interest),pch=19,
#         main=paste(file_name,' R^2=',sprintf('%0.3f',lm_sum$adj.r.squared),sep=''))
#     graphics.off()
# }
# 
# for (prop in data_of_interest) {
#     pdf(file.path(out_folder,paste('boxplots_',prop,'.pdf',sep='')))
#     boxplot(IA32_Pax$assem[[prop]],IA32_Vin$assem[[prop]],KD_Pax$assem[[prop]],KD_Vin$assem[[prop]],
#         names=c('IA32 Pax','IA32 Vin','KD Pax','KD Vin'),
#         notch=T,main=paste('Assembly',prop))
#     boxplot(IA32_Pax$dis[[prop]],IA32_Vin$dis[[prop]],KD_Pax$dis[[prop]],KD_Vin$dis[[prop]],
#         names=c('IA32 Pax','IA32 Vin','KD Pax','KD Vin'),
#         notch=T,main=paste('Disassembly',prop))
#     graphics.off()
# }

########################################
#Single Ad Plots
#######################################

# plot_all_ad_intensities(raw_data$IA32_Pax$intensity,IA32_Pax)
# plot_all_ad_intensities(raw_data$IA32_Vin$intensity,IA32_Vin)
# plot_all_ad_intensities(raw_data$KD_Pax$intensity,KD_Pax)
# plot_all_ad_intensities(raw_data$KD_Vin$intensity,KD_Vin)
 
###########################################################
#Kinetics Figures
###########################################################

#######################################
# Fibronectin Concentrations
#######################################

svg(file.path(out_folder,'Kinetics','fibro_concen_kinetics.svg'),width=8,height=5);
layout(cbind(1,2))
par(bty='n',mar=c(2.8,2.8,3,0), mgp=c(1.6,0.5,0),xpd=T)

#assembly plot
boxplot_with_points(list(Fibro_1ug$assembly$slope,Fibro_10ug$assembly$slope,Fibro_100ug$assembly$slope),
    names=c('1','10','100'),with.p.value=F,inc.points=F,
    ylab=expression(paste('FA Assembly Rate (',min^-1,')',sep='')))

pl_size = par("usr");

mtext('Fibronectin Concentration (\u03BCg/mL)',side=1,line=1.7)

plot_signif_bracket(c(1,pl_size[4]*0.93),c(2,pl_size[4]*0.92), 
    over_text=paste('p<', fibro_1ug_10ug_p$assembly$p.value[1],sep=''),text_sep=0)
plot_signif_bracket(c(2,pl_size[4]*0.99),c(3,pl_size[4]*0.98), 
    over_text=paste('p<',fibro_10ug_100ug_p$assembly$p.value[1],sep=''))
plot_signif_bracket(c(1,pl_size[4]*1.08),c(3,pl_size[4]*1.07), 
    over_text=paste('p<',fibro_1ug_100ug_p$assembly$p.value[1],sep=''))

#disassembly plot
boxplot_with_points(list(Fibro_1ug$dis$slope,Fibro_10ug$dis$slope,Fibro_100ug$dis$slope),
    names=c('1','10','100'),with.p.value=F,inc.points=F,
    ylab=expression(paste('FA Disassembly Rate (',min^-1,')',sep='')))

pl_size = par("usr");
mtext('Fibronectin Concentration (\u03BCg/mL)',side=1,line=1.7)
plot_signif_bracket(c(1,pl_size[4]*0.93),c(2,pl_size[4]*0.92), 
    over_text=paste('p<', fibro_1ug_10ug_p$dis$p.value[1],sep=''),text_sep=0)
plot_signif_bracket(c(2,pl_size[4]*0.99),c(3,pl_size[4]*0.98), 
    over_text=paste('p<',fibro_10ug_100ug_p$dis$p.value[1],sep=''))
plot_signif_bracket(c(1,pl_size[4]*1.08),c(3,pl_size[4]*1.07), 
    over_text=paste('p<',fibro_1ug_100ug_p$dis$p.value[1],sep=''))

graphics.off()

#######################################
# Different Adhesion Tags
#######################################

svg(file.path(out_folder,'Kinetics','all_tagged_kinetics.svg'),width=8,height=10);
layout(rbind(1,2))
par(bty='n',mar=c(2.8,2.8,2,0), mgp=c(1.6,0.5,0),xpd=T)

#assembly plot
boxplot_with_points(list(IA32_Pax$assembly$slope,KD_Pax$assembly$slope,
    IA32_FAK$assembly$slope,KD_FAK$assembly$slope,
    IA32_Vin$assembly$slope, KD_Vin$assembly$slope),
    with.p.value=F,inc.points=F, axes=F,
    ylab=expression(paste('FA Assembly Rate (',min^-1,')',sep='')))

pl_size = par("usr");
char_size = par("cxy")[2]
axis(2)

labs = c(paste('NS (n=',length(IA32_Pax$assembly$slope),')',sep=''),
    paste('2xKD (n=',length(KD_Pax$assembly$slope),')',sep=''))
axis(1, labels = labs, at = c(1,2))
lines(c(1,2),rep(pl_size[3]-char_size*1.675,2),lwd=3)
mtext('Pax',1,at=1.5,line=1.7)
plot_signif_bracket(c(1,pl_size[4]*1.01),c(2,pl_size[4]*1), 
    over_text=paste('p<', Pax_p$assembly$p.value[1],sep=''),text_sep=0)

labs = c(paste('NS (n=',length(IA32_FAK$assembly$slope),')',sep=''),
    paste('2xKD (n=',length(KD_FAK$assembly$slope),')',sep=''))
axis(1, labels = labs, at = c(3,4))
lines(c(3,4),rep(pl_size[3]-char_size*1.675,2),lwd=3)
mtext('FAK',1,at=3.5,line=1.7)
plot_signif_bracket(c(3,pl_size[4]*1.01),c(4,pl_size[4]*1), 
    over_text=paste('p<', FAK_p$assembly$p.value[1],sep=''),text_sep=0)

labs = c(paste('NS (n=',length(IA32_Vin$assembly$slope),')',sep=''),
    paste('2xKD (n=',length(KD_Vin$assembly$slope),')',sep=''))
axis(1, labels = labs, at = c(5,6))
lines(c(5,6),rep(pl_size[3]-char_size*1.675,2),lwd=3)
mtext('Vin',1,at=5.5,line=1.7)
plot_signif_bracket(c(5,pl_size[4]*1.01),c(6,pl_size[4]*1), 
    over_text=paste('p<', Vin_p$assembly$p.value[1],sep=''),text_sep=0)

#disassembly plot
boxplot_with_points(list(IA32_Pax$dis$slope,KD_Pax$dis$slope,
    IA32_FAK$dis$slope,KD_FAK$dis$slope,
    IA32_Vin$dis$slope, KD_Vin$dis$slope),
    with.p.value=F,inc.points=F, axes=F,
    ylab=expression(paste('FA Disassembly Rate (',min^-1,')',sep='')))

pl_size = par("usr");
char_size = par("cxy")[2]
axis(2)

labs = c(paste('NS (n=',length(IA32_Pax$dis$slope),')',sep=''),
    paste('2xKD (n=',length(KD_Pax$dis$slope),')',sep=''))
axis(1, labels = labs, at = c(1,2))
lines(c(1,2),rep(pl_size[3]-char_size*1.675,2),lwd=3)
mtext('Pax',1,at=1.5,line=1.7)
plot_signif_bracket(c(1,pl_size[4]*1.01),c(2,pl_size[4]*1), 
    over_text=paste('p<', Pax_p$dis$p.value[1],sep=''),text_sep=0)

labs = c(paste('NS (n=',length(IA32_FAK$dis$slope),')',sep=''),
    paste('2xKD (n=',length(KD_FAK$dis$slope),')',sep=''))
axis(1, labels = labs, at = c(3,4))
lines(c(3,4),rep(pl_size[3]-char_size*1.675,2),lwd=3)
mtext('FAK',1,at=3.5,line=1.7)
plot_signif_bracket(c(3,pl_size[4]*1.01),c(4,pl_size[4]*1), 
    over_text=paste('p<', FAK_p$dis$p.value[1],sep=''),text_sep=0)

labs = c(paste('NS (n=',length(IA32_Vin$dis$slope),')',sep=''),
    paste('2xKD (n=',length(KD_Vin$dis$slope),')',sep=''))
axis(1, labels = labs, at = c(5,6))
lines(c(5,6),rep(pl_size[3]-char_size*1.675,2),lwd=3)
mtext('Vin',1,at=5.5,line=1.7)
plot_signif_bracket(c(5,pl_size[4]*1.01),c(6,pl_size[4]*1), 
    over_text=paste('p<', Vin_p$dis$p.value[1],sep=''),text_sep=0)

graphics.off()

###########################################################
# FA Alignment Comparisons
###########################################################

source('FA_alignment_search.R')
svg(file.path(out_folder,'Kinetics','angle_split.svg'),width=7,height=3);
layout(cbind(1,2,3))
par(bty='n',mar=c(1.8,2.8,0,0), mgp=c(1.6,0.5,0),xpd=T)
boxplot_with_points(list(split_data$Fibro_100ug$middle_assembly$assembly_rate,
    split_data$Fibro_100ug$out_assembly$assembly_rate),notch=T,
    names=c('Middle','Out'),with.p.value=F, inc.points=F,
    p.value.type = 'median',ylab=expression(paste('FA Assembly Rate (',min^-1,')',sep='')))

boxplot_with_points(list(split_data$Fibro_100ug$middle_disassembly$disassembly_rate,
    split_data$Fibro_100ug$out_disassembly$disassembly_rate),notch=T,
    names=c('Middle','Out'),with.p.value=F,inc.points=F,
    p.value.type = 'median',ylab=expression(paste('FA Disassembly Rate (',min^-1,')',sep='')))

boxplot_with_points(list(split_data$Fibro_100ug$middle_longev$longevity*2.5,
    split_data$Fibro_100ug$out_longev$longevity*2.5),notch=T,
    names=c('Middle','Out'),with.p.value=F,inc.points = F,
    ylab='Longevity (min)')
graphics.off()
