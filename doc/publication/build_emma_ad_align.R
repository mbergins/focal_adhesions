rm(list = ls())
library(lattice)
library(geneplotter)
library(Hmisc)
source('FA_analysis_lib.R')
source('bilinear_modeling.R')
source('FA_alignment_search.R')

###############################################################################
# Results Gathering
###############################################################################

raw_data_base_dir = '../../results/emma/processed_2stdev/';

#######################################
# Fibronectin Concentrations
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_trial_2/WT_*/*/FA_orientation.Rdata*'))
fibro_1ug = load_alignment_props(align_files);

align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_trial_2/WT_*/*/FA_orientation.Rdata*'))
fibro_10ug = load_alignment_props(align_files);

align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_trial_2/WT_*/*/FA_orientation.Rdata*'))
fibro_100ug = load_alignment_props(align_files);

#######################################
# FAK WT/KD
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/FAK*/*/FA_orientation.Rdata*'))
IA32_FAK = load_alignment_props(align_files);

align_files = Sys.glob(file.path(raw_data_base_dir,'100ug_KD/FAK*/*/FA_orientation.Rdata*'))
KD_FAK = load_alignment_props(align_files);

#######################################
# Paxillin WT/KD
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/Pax*/*/FA_orientation.Rdata*'))
IA32_Pax = load_alignment_props(align_files);

align_files = Sys.glob(file.path(raw_data_base_dir,'100ug_KD/Pax*/*/FA_orientation.Rdata*'))
KD_Pax = load_alignment_props(align_files);

#######################################
# Vinculin WT/KD
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'100ug_IA32/Vin*/*/FA_orientation.Rdata*'))
IA32_Vin = load_alignment_props(align_files);

align_files = Sys.glob(file.path(raw_data_base_dir,'100ug_KD/Vin*/*/FA_orientation.Rdata*'))
KD_Vin = load_alignment_props(align_files);

#######################################
# Haptotaxis
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'Hapto/WT/Pax*/*/FA_orientation.Rdata'))
hapto_wt = load_alignment_props(align_files);
 
# align_files = Sys.glob(file.path(raw_data_base_dir,'Hapto/2xKD/Pax*/*/FA_orientation*'))
# hapto_kd = load_alignment_props(align_files)

system('notify-send "done with R"')
out_folder = '../../doc/publication/figures/emma/'

###########################################################
# P-value calculations
###########################################################

# fibro_1_10_p = determine_mean_p_value(fibro_1ug$best_FAAI,fibro_10ug$best_FAAI)
# fibro_1_100_p = determine_mean_p_value(fibro_1ug$best_FAAI,fibro_100ug$best_FAAI)
# fibro_10_100_p = determine_mean_p_value(fibro_10ug$best_FAAI,fibro_100ug$best_FAAI)
# 
# FAK_p = determine_mean_p_value(IA32_FAK$best_FAAI,KD_FAK$best_FAAI)
# Pax_p = determine_mean_p_value(IA32_Pax$best_FAAI,KD_Pax$best_FAAI)
# Vin_p = determine_mean_p_value(IA32_Vin$best_FAAI,KD_Vin$best_FAAI)

stop()
###############################################################################
# Plotting
###############################################################################

###########################################################
# FAAI
###########################################################

#######################################
# Paxillin in NS on Fibro
#######################################

svg(file.path(out_folder,'FAAI','fibro_FAAI_temp.svg'),width=4,height=4);

#the xpd=T part allows us to draw lines outside the plotting area, this is
#needed to get the lines under the lables
par(bty='n',mar=c(2.7,2.6,2.5,0), mgp=c(1.6,0.5,0),xpd=T)

boxplot_with_points(list(fibro_1ug$best_FAAI, fibro_10ug$best_FAAI, 
    fibro_100ug$best_FAAI,hapto_wt$best_FAAI,KD_Pax$best_FAAI),
    ylab='FA Alignment Index',with.p.value=F, names = c('1','10','100','Hapto','KD'),
    colors=c('black'))

boxplot_with_points(list(fibro_1ug$best_FAAI, fibro_10ug$best_FAAI, fibro_100ug$best_FAAI),
    ylab='FA Alignment Index',with.p.value=F, names = c('1','10','100'),
    colors=c('black'))

p_size = par("usr");
char_size = par("cxy")[2];

mtext('Fibronectin Concentration (\u03BCg/mL)',1,at=mean(p_size[1:2]),line=1.6)

#put in the significance brackets
plot_signif_bracket(c(1,p_size[4]*0.95),c(2,p_size[4]*0.945), 
    over_text=paste('p<',fibro_1_10_p$p.value[1],sep=''))
plot_signif_bracket(c(2,p_size[4]*1),c(3,p_size[4]*0.995), 
    over_text=paste('p<',fibro_10_100_p$p.value[1],sep=''))
plot_signif_bracket(c(1,p_size[4]*1.04),c(3,p_size[4]*1.035), 
    over_text=paste('p<',fibro_1_100_p$p.value[1],sep=''))
graphics.off()

#######################################
# Adhesion Tags in NS and 2xKD
#######################################

#placing the labels on this figure has become fairly complicated
svg(file.path(out_folder,'FAAI','tagged_ads_100ug_comparisons.svg'),width=7,height=4);

#the xpd=T part allows us to draw lines outside the plotting area, this is
#needed to get the lines under the lables
par(bty='n',mar=c(2.6,2.6,0,0), mgp=c(1.6,0.5,0),xpd=T)

#skip drawing axes
boxplot_with_points(list(IA32_Pax$best_FAAI,KD_Pax$best_FAAI,
    IA32_FAK$best_FAAI,KD_FAK$best_FAAI,
    IA32_Vin$best_FAAI,KD_Vin$best_FAAI),
    ylab='FA Alignment Index',ylim=c(40,80),with.p.value=F, 
    colors=c('black'),axes=F)

p_size = par("usr");
#draw the axes
axis(2)

labs = c(paste('NS (n=',length(IA32_Pax$best_FAAI),')',sep=''),
    paste('2xKD (n=',length(KD_Pax$best_FAAI),')',sep=''))
axis(1, labels = labs, at = c(1,2))
lines(c(1,2),rep(p_size[3]*0.89,2),lwd=3)
mtext('Pax',1,at=1.5,line=1.6)

labs = c(paste('NS (n=',length(IA32_FAK$best_FAAI),')',sep=''),
    paste('2xKD (n=',length(KD_FAK$best_FAAI),')',sep=''))
axis(1, labels = labs, at = c(3,4))
lines(c(3,4),rep(p_size[3]*0.89,2),lwd=3)
mtext('FAK',1,at=3.5,line=1.6)

labs = c(paste('NS (n=',length(IA32_Vin$best_FAAI),')',sep=''),
    paste('2xKD (n=',length(KD_Vin$best_FAAI),')',sep=''))
axis(1, labels = labs, at = c(5,6))
lines(c(5,6),rep(p_size[3]*0.89,2),lwd=3)
mtext('Vin',1,at=5.5,line=1.6)

#put in the significance brackets
plot_signif_bracket(c(1,p_size[4]*0.97),c(2,p_size[4]*0.965), 
    over_text=paste('p<',Pax_p$p.value[1],sep=''))
plot_signif_bracket(c(3,p_size[4]*0.97),c(4,p_size[4]*0.965), 
    over_text=paste('p<',FAK_p$p.value[1],sep=''))
plot_signif_bracket(c(5,p_size[4]*0.97),c(6,p_size[4]*0.965), 
    over_text=paste('p<',Vin_p$p.value[1],sep=''))
graphics.off()

#######################################
# Single Adhesion Stability
#######################################

source('FA_alignment_search.R')
align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_trial_2/WT_*/*/FA_orientation.Rdata'));
lin_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_trial_2/WT_*/adhesion_props/single_lin.csv'));
ad_dev_1ug = gather_all_single_adhesion_deviances(align_files,lin_files);

source('FA_alignment_search.R')
align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_trial_2/WT_*/*/FA_orientation.Rdata'));
lin_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_trial_2/WT_*/adhesion_props/single_lin.csv'));
ad_dev_10ug = gather_all_single_adhesion_deviances(align_files,lin_files);

source('FA_alignment_search.R')
align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_trial_2/WT_*/*/FA_orientation.Rdata'));
lin_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_trial_2/WT_*/adhesion_props/single_lin.csv'));
ad_dev_100ug = gather_all_single_adhesion_deviances(align_files,lin_files,min.area=-Inf, min.data.points=2);

svg(file.path(out_folder,'FAAI','single_ad_stability.svg'));
par(bty='n',mar=c(2.8,2.6,0.5,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
hist_100ug = hist(ad_dev_100ug$mean_dev,xlab='Average Deviance from Start (degrees)',
	ylab='# of Focal Adhesions',main='',axis=F)
axis(1,lwd=3)
axis(2,lwd=3)
graphics.off()

#######################################
# Simple FAAI figure hist
#######################################

fibro_100ug_10 = get(load('../../results/emma/processed_2stdev/Fibro_100ug_trial_2/WT_10/adhesion_props/FA_orientation.Rdata'));

svg(file.path(out_folder,'simple orientation','FAAI_hist.svg'),width=5,height=5);
par(bty='n',mar=c(2.7,2.5,0.75,0.1), mgp=c(1.5,0.5,0),xpd=T,cex=2,lwd=3)
hist(fibro_100ug_10$corr,axes=F,xlim=c(-90,90),ylab='Number of Adhesions',
    xlab='FA Orientation',main='');
axis(2,lwd=3);
axis(1,at=c(-90,0,90),lwd=3);
pl_size = par("usr");
FAAI = find_FAAI_from_orientation(fibro_100ug_10$corrected);
#text(pl_size[2],pl_size[4],'FAAI=90-stdev(FA Angles)',pos=2)
text(pl_size[2],pl_size[4],paste('FAAI=',sprintf('%.1f',FAAI),sep=''),pos=2)
graphics.off()

#######################################
# Sample FAAI plots
#######################################

var_name = load(as.character(fibro_100ug$align_file[1]))
high_FAAI_data = get(var_name)

svg(file.path(out_folder,'FAAI','sample_high_FAAI.svg'));
par(bty='n',mar=c(2.8,2.6,0.5,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
breaks = hist(high_FAAI_data$corrected_orientation,axes=F,ylab='Number of Adhesions',
    xlab='FA Orientation',main='')
axis(2, lwd=3)
axis(1,at=seq(-90,90,by=45),lwd=3)
pl_size = par("usr");
text(pl_size[2],pl_size[4],paste('FAAI=',sprintf('%.1f',fibro_100ug$best_FAAI[1]),sep=''),pos=2)
graphics.off()

var_name = load(as.character(fibro_1ug$align_file[8]))
low_FAAI_data = get(var_name)

svg(file.path(out_folder,'FAAI','sample_low_FAAI.svg'));
par(bty='n',mar=c(2.6,2.8,0.5,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
breaks = hist(low_FAAI_data$corrected_orientation,axes=F,ylab='Number of Adhesions',
    xlab='FA Orientation',main='')
axis(2, lwd=3)
axis(1,at=seq(-90,90,by=45),lwd=3)
pl_size = par("usr");
text(pl_size[2],pl_size[4],paste('FAAI=',sprintf('%.1f',fibro_1ug$best_FAAI[8]),sep=''),pos=2)
graphics.off()

#######################################
# Sample FA orientation cartoon
#######################################

# load("../../results/emma/Fibro_100ug/WT_02/adhesion_props/FA_orientation.Rdata"); 
# angle_search = test_dom_angles(data_set$subseted_data$orientation);
# best_angle = find_best_alignment_angle(angle_search)
# corrected_orientation = apply_new_orientation(data_set$subseted_data$orientation,best_angle)
# 
# per_image_dom_angle = find_per_image_dom_angle(data_set$mat, min_eccen=3)
# 
# deg_45_orientation = apply_new_orientation(data_set$subseted_data$orientation,45)
# deg_90_orientation = apply_new_orientation(data_set$subseted_data$orientation,90)
# deg_135_orientation = apply_new_orientation(data_set$subseted_data$orientation,135)
# deg_179_orientation = apply_new_orientation(data_set$subseted_data$orientation,179)
# 
# or_data_cartoon = list(data = rbind(data_set$subseted_data$orientation,
#     deg_45_orientation,deg_90_orientation,deg_135_orientation),
#     filenames = c('deg_0_hist.svg','deg_45_hist.svg','deg_90_hist.svg','deg_135_hist.svg'));
# 
# for (i in 1:length(or_data_cartoon$filenames)) {
#     svg(file.path(out_folder,'orientation_cartoon',or_data_cartoon$filenames[i]),pointsize=36)
#     par(bty='n',mar=c(2.5,2.25,0.0,0.5), mgp=c(1.4,0.5,0))
#     hist(or_data_cartoon$data[i,],breaks=seq(-90,90,by=10),xlim=c(-90,90),axes=F,main='',
#         xlab='Adhesion Angle',col='black')
#     axis(2,lwd=3)
#     axis(1,at=c(-90,0,90),lwd=3)
#     if (i < 3) {
#         text(-90,550,paste('FAAI=',sprintf('%.1f',90-sd(or_data_cartoon$data[i,]))),col='red',pos=4)
#     } else {
#         text(90,550,paste('FAAI=',sprintf('%.1f',90-sd(or_data_cartoon$data[i,]))),col='red',pos=2)
#     }
#     graphics.off()
# }
# 
# svg(file.path(out_folder,'orientation_cartoon','full_FAAI_plot.svg'),height=4,width=14,pointsize=24)
# 
# par(bty='n',mar=c(2.5,2.25,0.8,0.5), mgp=c(1.4,0.5,0))
# plot(angle_search$test_angles,angle_search$angle_FAAI,typ='l',lwd=3,xlab='Alignment Angles (degrees)',ylab='Alignment Index',axes=F)
# axis(2,lwd=3)
# axis(1,lwd=3,at=seq(0,180,by=20))
# 
# highlight_points = c(0,45,90,135);
# for (angle in highlight_points) {
#     degree_index = which(angle_search$test_angles == angle)
#     points(angle_search$test_angles[degree_index],
#         angle_search$angle_FAAI[degree_index],col='red',pch=19)
# }
# 
# points(best_angle, 90-sd(corrected_orientation),col='green',pch=19)
# 
# graphics.off()
