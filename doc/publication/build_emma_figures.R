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
# Various adhesion tags, all on 100ug fibronectin
########################################

# #KD - Pax
# exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_KD/Pax_*/adhesion_props/models/'))
# exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
# raw_data$KD_Pax$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));
# 
# #IA32 - Pax
# exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/Pax_*/adhesion_props/models/'))
# exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
# raw_data$WT_Pax$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));
# 
# #KD - FAK 
# exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_KD/FAK_*/adhesion_props/models/'))
# exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
# raw_data$KD_FAK$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));
# 
# #IA32 - FAK 
# exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/FAK_*/adhesion_props/models/'))
# exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
# raw_data$WT_FAK$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));
# 
# #KD - Vin
# exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_KD/Vin_*/adhesion_props/models/'))
# exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
# raw_data$KD_Vin$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));
# 
# #IA32 - Vin
# exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/Vin_*/adhesion_props/models/'))
# exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
# raw_data$WT_Vin$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));

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

#######################################
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
