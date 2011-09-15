rm(list = ls())
library(lattice)
library(geneplotter)
library(Hmisc)
source('FA_analysis_lib.R')
source('bilinear_modeling.R')
source('FA_alignment_search.R')
source('../../doc/publication/errbar.s')

###############################################################################
# Results Gathering
###############################################################################

raw_data_base_dir = '../../results/emma/processed_2stdev/';

#######################################
# Fibronectin Concentrations - 2xKD
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_2xKD/Pax_*/*/FA_orientation.Rdata*'))
fibro_1ug_2xKD = load_alignment_props(align_files);

align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_2xKD/Pax_*/*/FA_orientation.Rdata*'))
fibro_10ug_2xKD = load_alignment_props(align_files);

align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_2xKD/Pax_*/*/FA_orientation.Rdata*'))
fibro_100ug_2xKD = load_alignment_props(align_files);

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
# Rat 2 Control
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'Rat2/*/*/FA_orientation.Rdata*'))
rat2 = load_alignment_props(align_files);

#######################################
# Rat 2 - CK666
#######################################

align_files = Sys.glob(file.path(raw_data_base_dir,'Rat2_ck666/*/*/FA_orientation.Rdata*'))
rat2_ck666 = load_alignment_props(align_files);

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

out_folder = '../../doc/publication/figures/emma/'

###########################################################
# P-value calculations
###########################################################

fibro_1_10_p = determine_mean_p_value(fibro_1ug$best_FAAI,fibro_10ug$best_FAAI)
fibro_1_100_p = determine_mean_p_value(fibro_1ug$best_FAAI,fibro_100ug$best_FAAI)
fibro_10_100_p = determine_mean_p_value(fibro_10ug$best_FAAI,fibro_100ug$best_FAAI)

fibro_1_10_KD_p = determine_mean_p_value(fibro_1ug_2xKD$best_FAAI,fibro_10ug_2xKD$best_FAAI)
fibro_1_100_KD_p = determine_mean_p_value(fibro_1ug_2xKD$best_FAAI,fibro_100ug_2xKD$best_FAAI)
fibro_10_100_KD_p = determine_mean_p_value(fibro_10ug_2xKD$best_FAAI,fibro_100ug_2xKD$best_FAAI)

fibro_NS_KD_1ug_p = determine_mean_p_value(fibro_1ug_2xKD$best_FAAI,fibro_1ug$best_FAAI)
fibro_NS_KD_10ug_p = determine_mean_p_value(fibro_10ug_2xKD$best_FAAI,fibro_10ug$best_FAAI)
fibro_NS_KD_100ug_p = determine_mean_p_value(fibro_100ug_2xKD$best_FAAI,fibro_100ug$best_FAAI)

rat2_p = determine_mean_p_value(rat2$best_FAAI,rat2_ck666$best_FAAI)

Pax_p = determine_mean_p_value(IA32_Pax$best_FAAI,KD_Pax$best_FAAI)
FAK_p = determine_mean_p_value(IA32_FAK$best_FAAI,KD_FAK$best_FAAI)
Vin_p = determine_mean_p_value(IA32_Vin$best_FAAI,KD_Vin$best_FAAI)

large_area_p = determine_mean_p_value(fibro_100ug_2xKD$large_FAAI,fibro_100ug$large_FAAI)
medium_area_p = determine_mean_p_value(fibro_100ug_2xKD$medium_FAAI,fibro_100ug$medium_FAAI)
small_area_p = determine_mean_p_value(fibro_100ug_2xKD$small_FAAI,fibro_100ug$small_FAAI)

small_medium_NS = determine_mean_p_value(fibro_100ug$small_FAAI,fibro_100ug$medium_FAAI)
small_large_NS = determine_mean_p_value(fibro_100ug$small_FAAI,fibro_100ug$large_FAAI)
medium_large_NS = determine_mean_p_value(fibro_100ug$medium_FAAI,fibro_100ug$large_FAAI)

small_medium_2xKD = determine_mean_p_value(fibro_100ug_2xKD$small_FAAI,fibro_100ug_2xKD$medium_FAAI)
small_large_2xKD = determine_mean_p_value(fibro_100ug_2xKD$small_FAAI,fibro_100ug_2xKD$large_FAAI)
medium_large_2xKD = determine_mean_p_value(fibro_100ug_2xKD$medium_FAAI,fibro_100ug_2xKD$large_FAAI)

system('notify-send "done with R"')
stop()
###############################################################################
# Plotting
###############################################################################

KD_2X_color = rgb(148/255,150/255,152/255,1);

###########################################################
# FAAI
###########################################################

#######################################
# Paxillin in NS on Fibro
#######################################

svg(file.path(out_folder,'FAAI','fibro_FAAI.svg'),width=3.1,height=4);

#the xpd=T part allows us to draw lines outside the plotting area, this is
#needed to get the lines under the lables
par(bty='n',mar=c(2.5,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T,lwd=1.5)

boxplot(fibro_1ug$best_FAAI,fibro_1ug_2xKD$best_FAAI, 
    fibro_10ug$best_FAAI,fibro_10ug_2xKD$best_FAAI, 
    fibro_100ug$best_FAAI,fibro_100ug_2xKD$best_FAAI,
    ylab='FA Alignment Index',axes=F,col=c('white',KD_2X_color))

axis(2,lwd=1.5)
axis(1,at=c(1.5,3.5,5.5),labels=c(1,10,100),lwd=1.5)

legend('topleft',c('NS','2xKD'),fill=c('white',KD_2X_color),bty='n')

p_size = par("usr")

plot_signif_bracket(c(1,p_size[4]*0.875),c(2,p_size[4]*0.87), 
    over_text=paste('p<',fibro_NS_KD_1ug_p$p.value[1],sep=''))

plot_signif_bracket(c(3,p_size[4]*0.94),c(4,p_size[4]*0.935), 
    over_text=paste('p<',fibro_NS_KD_10ug_p$p.value[1],sep=''))

plot_signif_bracket(c(5,p_size[4]*1),c(6,p_size[4]*0.995), 
    over_text=paste('p<',fibro_NS_KD_100ug_p$p.value[1],sep=''))

mtext('Fibronectin Concentration (\u03BCg/mL)',1,at=mean(p_size[1:2]),line=1.4)

graphics.off()

#######################################
# Adhesion Tags in NS and 2xKD
#######################################

#placing the labels on this figure has become fairly complicated
svg(file.path(out_folder,'FAAI','tagged_ads_100ug_comparisons.svg'),width=4.75,height=4);

#the xpd=T part allows us to draw lines outside the plotting area, this is
#needed to get the lines under the lables
par(bty='n',mar=c(3.7,2.6,1,0), mgp=c(1.6,0.5,0),xpd=T,lwd=1.5)

#skip drawing axes
boxplot(IA32_Pax$best_FAAI,KD_Pax$best_FAAI,
    IA32_FAK$best_FAAI,KD_FAK$best_FAAI,
    IA32_Vin$best_FAAI,KD_Vin$best_FAAI,
    ylab='FA Alignment Index',axes=F,col=c('white',KD_2X_color))

p_size = par("usr");
#draw the axes
axis(2,lwd=1.5)

# legend('topleft',c('NS','2xKD'),fill=c('white',KD_2X_color),bty='n')

text_height = strheight('Paxillin')*0.25
mtext('Paxillin',at=c(1.5,p_size[4]))
segments(1,p_size[4]-text_height,2,p_size[4]-text_height,lwd=3)
mtext('FAK',at=c(3.5,p_size[4]))
segments(3,p_size[4]-text_height,4,p_size[4]-text_height,lwd=3)
mtext('Vinculin',at=c(5.5,p_size[4]))
segments(5,p_size[4]-text_height,6,p_size[4]-text_height,lwd=3)

labs = c(paste('NS \n(n=',length(IA32_Pax$best_FAAI),')',sep=''),
    paste('2xKD \n(n=',length(KD_Pax$best_FAAI),')',sep=''))
axis(1, labels = labs, at = c(1,2),padj=0.6,lwd=1.5)
mtext(paste('p<',Pax_p$p.value[1],sep=''),side=1,at=c(1.5,p_size[3]),line=2.6)

labs = c(paste('NS \n(n=',length(IA32_FAK$best_FAAI),')',sep=''),
    paste('2xKD \n(n=',length(KD_FAK$best_FAAI),')',sep=''))
axis(1, labels = labs, at = c(3,4),padj=0.6,lwd=1.5)
mtext(paste('p<',FAK_p$p.value[1],sep=''),side=1,at=c(3.5,p_size[3]),line=2.6)
 
labs = c(paste('NS \n(n=',length(IA32_Vin$best_FAAI),')',sep=''),
    paste('2xKD \n(n=',length(KD_Vin$best_FAAI),')',sep=''))
axis(1, labels = labs, at = c(5,6),padj=0.6,lwd=1.5)
mtext(paste('p<',Vin_p$p.value[1],sep=''),side=1,at=c(5.5,p_size[3]),line=2.6)

graphics.off()

#######################################
# Paxillin in Rat2
#######################################
svg(file.path(out_folder,'FAAI','Rat2_FAAI.svg'),width=2.2,height=4);

#the xpd=T part allows us to draw lines outside the plotting area, this is
#needed to get the lines under the lables
par(bty='n',mar=c(2.6,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T,lwd=1.5)

rat2_names = c(paste('CK-689 \n(n=',dim(rat2)[1],')',sep=''),paste('CK-666 \n(n=',dim(rat2_ck666)[1],')',sep=''))

boxplot(rat2$best_FAAI,rat2_ck666$best_FAAI,
    ylab='FA Alignment Index',axes=F)

axis(2,lwd=1.5);
axis(1,labels = rat2_names,at=c(1,2),padj=0.6,lwd=1.5);

p_size = par("usr");

plot_signif_bracket(c(1,p_size[4]*1),c(2,p_size[4]*0.995), 
    over_text=paste('p<',rat2_p$p.value[1],sep=''))

graphics.off()

#######################################
# Area Bins (NS vs 2xKD)
#######################################
svg(file.path(out_folder,'FAAI','area_bins_FAAI.svg'),width=4.75,height=4);

#the xpd=T part allows us to draw lines outside the plotting area, this is
#needed to get the lines under the lables
par(bty='n',mar=c(2.6,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T,lwd=1.5)

boxplot(fibro_100ug$small_FAAI,fibro_100ug_2xKD$small_FAAI,
    fibro_100ug$medium_FAAI,fibro_100ug_2xKD$medium_FAAI,
    fibro_100ug$large_FAAI,fibro_100ug_2xKD$large_FAAI,
    axes=F,ylab='FA Alignment Index',col=c('white',KD_2X_color))

axis(2)

p_size = par("usr");

text_height = strheight('Paxillin')*0.25
mtext('Small',at=c(1.5,p_size[4]))
segments(1,p_size[4]-text_height,2,p_size[4]-text_height,lwd=3)
mtext('Medium',at=c(3.5,p_size[4]))
segments(3,p_size[4]-text_height,4,p_size[4]-text_height,lwd=3)
mtext('Large',at=c(5.5,p_size[4]))
segments(5,p_size[4]-text_height,6,p_size[4]-text_height,lwd=3)

labs = c(paste('NS \n(n=',length(fibro_100ug$small_FAAI),')',sep=''),
    paste('2xKD \n(n=',length(fibro_100ug_2xKD$small_FAAI),')',sep=''))
axis(1, labels = labs, at = c(1,2),padj=0.6,lwd=1.5)
# mtext(paste('p<',Pax_p$p.value[1],sep=''),side=1,at=c(1.5,p_size[3]),line=2.6)

labs = c(paste('NS \n(n=',length(fibro_100ug$medium_FAAI),')',sep=''),
    paste('2xKD \n(n=',length(fibro_100ug_2xKD$medium_FAAI),')',sep=''))
axis(1, labels = labs, at = c(3,4),padj=0.6,lwd=1.5)
# mtext(paste('p<',FAK_p$p.value[1],sep=''),side=1,at=c(3.5,p_size[3]),line=2.6)
 
labs = c(paste('NS \n(n=',length(fibro_100ug$large_FAAI),')',sep=''),
    paste('2xKD \n(n=',length(fibro_100ug_2xKD$large_FAAI),')',sep=''))
axis(1, labels = labs, at = c(5,6),padj=0.6,lwd=1.5)
# mtext(paste('p<',Vin_p$p.value[1],sep=''),side=1,at=c(5.5,p_size[3]),line=2.6)

graphics.off()

###########################################################
# Simple FAAI Figure hist
###########################################################
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

###########################################################
# Sample FAAI plots
###########################################################
high_FAAI_data = get(load('../../results/emma/processed_2stdev/Fibro_100ug_trial_2/WT_01/adhesion_props/FA_orientation.Rdata'));

svg(file.path(out_folder,'sample_FAAI_stills','sample_high_FAAI.svg'));
par(bty='n',mar=c(2.6,2.6,0.6,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
breaks = hist(high_FAAI_data$corrected_orientation,axes=F,ylab='Number of Adhesions',
    xlab='FA Orientation',main='',xlim=c(-90,90),breaks=seq(-90,90,by=15))
axis(2, lwd=3)
axis(1,at=seq(-90,90,by=45),lwd=3)
pl_size = par("usr");
text(90,pl_size[4]*0.9,paste('FAAI=',sprintf('%.1f',fibro_100ug$best_FAAI[1]),sep=''),pos=2,cex=1.5)
graphics.off()


low_FAAI_data = get(load('../../results/emma/processed_2stdev//Fibro_1ug_trial_2/WT_08/adhesion_props/FA_orientation.Rdata'));

svg(file.path(out_folder,'sample_FAAI_stills','sample_low_FAAI.svg'));
par(bty='n',mar=c(2.6,2.6,0,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
hist(low_FAAI_data$corrected_orientation,axes=F,ylab='Number of Adhesions',
    xlab='FA Orientation',main='',xlim=c(-90,90),breaks=seq(-90,90,by=20))
axis(2, lwd=3)
axis(1,at=seq(-90,90,by=45),lwd=3)
pl_size = par("usr");
text(-90,pl_size[4]*0.9,paste('FAAI=',sprintf('%.1f',fibro_1ug$best_FAAI[8]),sep=''),pos=4,cex=1.5)
graphics.off()

###########################################################
# Single Adhesion Stability
###########################################################
align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_trial_2/WT_*/*/FA_orientation.Rdata'));
ad_dev_100ug = list();
for (file in align_files) {
    temp = get(load(file))
    temp = temp$single_ad_deviances
    temp$mean_area = temp$lineage$mean_area
    temp$longevity = temp$lineage$longevity
    temp$mean_axial_ratio = temp$lineage$mean_axial_ratio
    temp$exp_file = file
    temp = subset(temp, select=-c(lineage))
    ad_dev_100ug = rbind(ad_dev_100ug,temp)
}

align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_2xKD/Pax*/*/FA_orientation.Rdata*'))
ad_dev_100ug_2xKD = list();
for (file in align_files) {
    temp = get(load(file))
    temp = temp$single_ad_deviances
    temp$mean_area = temp$lineage$mean_area
    temp$longevity = temp$lineage$longevity
    temp$mean_axial_ratio = temp$lineage$mean_axial_ratio
    temp = subset(temp, select=-c(lineage))
    ad_dev_100ug_2xKD = rbind(ad_dev_100ug_2xKD,temp)
}

#############################
# Beside Histogram
#############################
hist_breaks = seq(0,90,by=5);

NS_hist_data = hist(ad_dev_100ug$mean_dev,plot=F,breaks=hist_breaks);
KD_hist_data = hist(ad_dev_100ug_2xKD$mean_dev,plot=F,breaks=hist_breaks);

all_hist_data = rbind(NS_hist_data$counts/(sum(NS_hist_data$counts)),
    KD_hist_data$counts/(sum(KD_hist_data$counts)))

bin_labels = c()
for (i in 1:length(hist_breaks[-1])) {
    bin_labels[i] = sprintf('%d - %d',hist_breaks[i],hist_breaks[i+1]);
}

svg(file.path(out_folder,'FAAI','single_ad_stability_dual.svg'), width=6,height=4);
par(bty='n',mar=c(4.4,2.5,0.5,0), mgp=c(1.6,0.5,0),xpd=T,lwd=1)
bar_data = barplot(all_hist_data,beside=T,axes=F,ylab='Percentage of Focal Adhesions',
    col=c('white',KD_2X_color),space=c(0,0.25))
axis(1,at = (bar_data[1,]+bar_data[2,])/2, labels=bin_labels,tick=F,las=2)
axis(2)
mtext('Average Deviance from Start (degrees)',side=1,line=3.25)

legend('bottomright',c('NS','2xKD'),
    fill=c('white',KD_2X_color),bty='n',inset=c(0.0,0.02),cex=1.5)

graphics.off()

#############################
# Overlapping Histogram
#############################
svg(file.path(out_folder,'FAAI','single_ad_stability_dual.svg'), width=10);
par(bty='n',mar=c(2.8,2.6,0.5,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
hist_data = hist(ad_dev_100ug$mean_dev,xlab='Average Deviance from Start (degrees)',
	ylab='Percentage of Focal Adhesions',main='',axes=F,lwd=3,col=rgb(0,1,0,0.5),
    freq=F,xlim=c(0,90))
hist(ad_dev_100ug_2xKD$mean_dev,axes=F,lwd=3,col=rgb(1,0,0,0.5),add=T,freq=F,breaks=seq(0,90,by=5))
axis(2,lwd=3)
axis(1,lwd=3,at=seq(0,90,by=30))

legend('bottomright',c('NS','2xKD'),cex=2,
    fill=c(rgb(0,1,0,0.5),rgb(1,0,0,0.5)),bty='n',inset=c(0.01,0.01))

# legend_text = c(paste('NS (n=',length(ad_dev_100ug$mean_dev),sep='',')')
#     ,paste('2xKD (n=',length(ad_dev_100ug_2xKD$mean_dev),sep='',')'))
# legend('bottomright',legend_text,
#     fill=c(rgb(0,1,0,0.2),rgb(1,0,0,0.2)),bty='n')

graphics.off()

###########################################################
# Spatial
###########################################################

#######################################
# Fibronectin Concentrations
#######################################
spatial_files = Sys.glob(file.path(raw_data_base_dir,'Fibro*trial_2/*/*/FA_dist_orientation.Rdata'))
all_spatial = list();

for (file in spatial_files) {
    this_data = get(load(file));
    if (regexpr('100ug',file) != -1) {
        print(paste(file,'100ug'));
        this_data$concen = rep('ug_100',dim(this_data)[1]);
    } else if (regexpr('10ug',file) != -1) {
        print(paste(file,'10ug'));
        this_data$concen = rep('ug_10',dim(this_data)[1]);
    } else if (regexpr('1ug',file) != -1) {
        print(paste(file,'1ug'));
        this_data$concen = rep('ug_1',dim(this_data)[1]);
    } else {
        print(paste('Couldnt find appropriate concentration',file))
    }
    all_spatial = rbind(all_spatial,this_data);
}

all_spatial$dist_bin = floor(all_spatial$dists/(2.5/0.1333));
concen_summary = list();
source('FA_alignment_search.R');
for (this_concen in unique(all_spatial$concen)) {
    this_concen_set = subset(all_spatial, concen==this_concen);
    concen_summary[[this_concen]]=bin_distance_data(this_concen_set);
}

svg(file.path(out_folder,'spatial','fibro_dist_vs_orientation.svg'),width=8,height=4);
layout(cbind(1,2))
par(bty='n',mar=c(2.7,2.6,0.6,0.2), mgp=c(1.6,0.5,0),xpd=T)
errbar(concen_summary$ug_1$dist, concen_summary$ug_1$or_mean,
    concen_summary$ug_1$or_plus,concen_summary$ug_1$or_minus,
    xlab='Distance Between Adhesions (\u03BCm)',ylab='Mean Angle Difference',
    col='green',xlim=c(0,80),ylim=c(0,45))
lines(lowess(concen_summary$ug_1$dist,concen_summary$ug_1$or_mean),col='green',lwd=3)

errbar(concen_summary$ug_100$dist, concen_summary$ug_100$or_mean,
    concen_summary$ug_100$or_plus,concen_summary$ug_100$or_minus,
    add=T)
lines(lowess(concen_summary$ug_100$dist,concen_summary$ug_100$or_mean),lwd=3)

errbar(concen_summary$ug_10$dist, concen_summary$ug_10$or_mean,
    concen_summary$ug_10$or_plus,concen_summary$ug_10$or_minus,
    add=T,col='red')
lines(lowess(concen_summary$ug_10$dist,concen_summary$ug_10$or_mean),col='red',lwd=3)

legend('bottomright',c('100ug','10ug','1ug'),fill=c('black','red','green'),inset=.01)

plot(concen_summary$ug_100$dist, concen_summary$ug_100$n_count,
    ylim=c(0,max(concen_summary$ug_100$n_count)),xlab='Distance Between Adhesions',
    ylab='Adhesion Count', typ='b')
lines(concen_summary$ug_10$dist, concen_summary$ug_10$n_count,col='red',typ='b')
lines(concen_summary$ug_1$dist, concen_summary$ug_1$n_count,col='green',typ='b')
graphics.off()

#######################################
# Different Adhesion Tags
#######################################
spatial_files = Sys.glob(file.path(raw_data_base_dir,'100ug*/*/*/FA_dist_orientation.Rdata'))
tag_spatial = list();

for (file in spatial_files) {
    this_data = get(load(file));
    
    output_str = c(file)
    if (regexpr('Pax',file) != -1) {
        this_data$tag = rep('Pax',dim(this_data)[1]);
        output_str = paste(output_str,'Pax');
    } else if (regexpr('Vin',file) != -1) {
        this_data$tag = rep('Vin',dim(this_data)[1]);
        output_str = paste(output_str,'Vin');
    } else if (regexpr('FAK',file) != -1) {
        this_data$tag = rep('FAK',dim(this_data)[1]);
        output_str = paste(output_str,'FAK');
    } else {
        print(paste('Cant find tag type',file))
    }
    
    if (regexpr('KD',file) != -1) {
        this_data$cell_type = rep('KD',dim(this_data)[1]);
        output_str = paste(output_str,'KD');
    } else if (regexpr('IA32',file) != -1) {
        this_data$cell_type = rep('IA32',dim(this_data)[1]);
        output_str = paste(output_str,'IA32');
    } else {
        print(paste('Couldnt find appropriate cell type',file))
    }
    tag_spatial = rbind(tag_spatial,this_data);
    print(output_str)
}

tag_spatial$dist_bin = floor(tag_spatial$dists/(2.5/0.1333));
tag_summary = list();
for (this_tag in unique(tag_spatial$tag)) {
    for (this_cell_type in unique(tag_spatial$cell_type)) {
        this_dist_data = subset(tag_spatial,tag == this_tag & cell_type == this_cell_type);
        tag_summary[[this_tag]][[this_cell_type]]=bin_distance_data(this_dist_data);
    }
}

svg(file.path(out_folder,'spatial','tags_dist_vs_orientation.svg'),width=8,height=12);
layout(rbind(c(1,2),c(3,4),c(5,6)))
for (tag_type in c('Pax','Vin','FAK')) {
    par(bty='n',mar=c(2.7,2.6,3,0.2), mgp=c(1.6,0.5,0),xpd=F)
    errbar(tag_summary[[tag_type]]$KD$dist, tag_summary[[tag_type]]$KD$or_mean,
        tag_summary[[tag_type]]$KD$or_plus,tag_summary[[tag_type]]$KD$or_minus,
        col='red',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',
        ylim=c(0,max(tag_summary[[tag_type]]$KD$or_plus)))
    lines(lowess(tag_summary[[tag_type]]$KD$dist,tag_summary[[tag_type]]$KD$or_mean),col='red',lwd=3)
    title(main=tag_type)

    errbar(tag_summary[[tag_type]]$IA32$dist, tag_summary[[tag_type]]$IA32$or_mean,
        tag_summary[[tag_type]]$IA32$or_plus,tag_summary[[tag_type]]$IA32$or_minus,
        col='green',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',add=T)
    lines(lowess(tag_summary[[tag_type]]$IA32$dist,tag_summary[[tag_type]]$IA32$or_mean),col='green',lwd=3)

    plot(tag_summary[[tag_type]]$IA32$dist, tag_summary[[tag_type]]$IA32$n_count,
        ylim=c(0,max(tag_summary[[tag_type]]$IA32$n_count)),xlab='Distance Between Adhesions',
        ylab='Adhesion Count', col='green',main=tag_type)
    lines(tag_summary[[tag_type]]$IA32$dist, tag_summary[[tag_type]]$IA32$n_count, col='green')
    points(tag_summary[[tag_type]]$KD$dist, tag_summary[[tag_type]]$KD$n_count,col='red')
    lines(tag_summary[[tag_type]]$KD$dist, tag_summary[[tag_type]]$KD$n_count,col='red')
}
graphics.off()

svg(file.path(out_folder,'spatial','tags_dist_vs_orientation_all.svg'));
errbar(tag_summary$Vin$KD$dist, tag_summary$Vin$KD$or_mean,
    tag_summary$Vin$KD$or_plus,tag_summary$Vin$KD$or_minus,
    col='red',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',
    ylim=c(0,max(tag_summary$Vin$KD$or_plus)),pch=1)
lines(lowess(tag_summary$Vin$KD$dist,tag_summary$Vin$KD$or_mean),col='red',lwd=3,lty=3)

errbar(tag_summary$Vin$IA32$dist, tag_summary$Vin$IA32$or_mean,
    tag_summary$Vin$IA32$or_plus,tag_summary$Vin$IA32$or_minus,
    col='red',add=T)
lines(lowess(tag_summary$Vin$IA32$dist,tag_summary$Vin$IA32$or_mean),col='red',lwd=3)

errbar(tag_summary$FAK$KD$dist, tag_summary$FAK$KD$or_mean,
    tag_summary$FAK$KD$or_plus,tag_summary$FAK$KD$or_minus,
    col='green',add=T,pch=1)
lines(lowess(tag_summary$FAK$KD$dist,tag_summary$FAK$KD$or_mean),col='green',lwd=3,lty=3)

errbar(tag_summary$FAK$IA32$dist, tag_summary$FAK$IA32$or_mean,
    tag_summary$FAK$IA32$or_plus,tag_summary$FAK$IA32$or_minus,
    col='green',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',add=T)
lines(lowess(tag_summary$FAK$IA32$dist,tag_summary$FAK$IA32$or_mean),col='green',lwd=3)

errbar(tag_summary$Pax$KD$dist, tag_summary$Pax$KD$or_mean,
    tag_summary$Pax$KD$or_plus,tag_summary$Pax$KD$or_minus,
    ,add=T,pch=1) 
lines(lowess(tag_summary$Pax$KD$dist,tag_summary$Pax$KD$or_mean),lwd=3,lty=3)

errbar(tag_summary$Pax$IA32$dist, tag_summary$Pax$IA32$or_mean,
    tag_summary$Pax$IA32$or_plus,tag_summary$Pax$IA32$or_minus,
    xlab='Distance Between Adhesions',ylab='Mean Angle Difference',add=T)
lines(lowess(tag_summary$Pax$IA32$dist,tag_summary$Pax$IA32$or_mean),lwd=3)

legend('bottomright',c('Pax','FAK','Vin'),fill=c('black','green','red'),inset=.01)
graphics.off()

###########################################################
# Sample FA orientation cartoon
###########################################################

cartoon_data = get(load("../../results/emma/processed_2stdev/Fibro_100ug/WT_02/adhesion_props/FA_orientation.Rdata")); 
angle_search = test_dom_angles(cartoon_data$subseted_data$orientation);
best_angle = find_best_alignment_angle(angle_search)
corrected_orientation = apply_new_orientation(cartoon_data$subseted_data$orientation,best_angle)

per_image_dom_angle = find_per_image_dom_angle(cartoon_data$mat, min_eccen=3)

deg_45_orientation = apply_new_orientation(cartoon_data$subseted_data$orientation,45)
deg_90_orientation = apply_new_orientation(cartoon_data$subseted_data$orientation,90)
deg_135_orientation = apply_new_orientation(cartoon_data$subseted_data$orientation,135)
deg_179_orientation = apply_new_orientation(cartoon_data$subseted_data$orientation,179)

or_data_cartoon = list(data = rbind(cartoon_data$subseted_data$orientation,
    deg_45_orientation,deg_90_orientation,deg_135_orientation),
    filenames = c('deg_0_hist.svg','deg_45_hist.svg','deg_90_hist.svg','deg_135_hist.svg'));

for (i in 1:length(or_data_cartoon$filenames)) {
    svg(file.path(out_folder,'orientation_cartoon',or_data_cartoon$filenames[i]),pointsize=36)
    par(bty='n',mar=c(2.5,2.25,0.0,0.5), mgp=c(1.4,0.5,0))
    hist(or_data_cartoon$data[i,],breaks=seq(-90,90,by=10),xlim=c(-90,90),axes=F,main='',
        xlab='Adhesion Angle',col='black')
    axis(2,lwd=3)
    axis(1,at=c(-90,0,90),lwd=3)
    if (i == 1) {
        text(-90,550,'FAAI=90-SD',col='red',pos=4)
        text(-90,490,paste('FAAI=',sprintf('%.1f',90-sd(or_data_cartoon$data[i,]))),col='red',pos=4)
    } else if (i == 2 || i == 3) {
        text(-90,550,paste('FAAI=',sprintf('%.1f',90-sd(or_data_cartoon$data[i,]))),col='red',pos=4)
    } else {
        text(90,550,paste('FAAI=',sprintf('%.1f',90-sd(or_data_cartoon$data[i,]))),col='red',pos=2)
    }
    graphics.off()
}

svg(file.path(out_folder,'orientation_cartoon','full_FAAI_plot.svg'),height=4,width=14,pointsize=24)

par(bty='n',mar=c(2.5,2.25,0.8,0.5), mgp=c(1.4,0.5,0))
plot(angle_search$test_angles,angle_search$angle_FAAI,typ='l',lwd=3,xlab='Alignment Angles (degrees)',ylab='Alignment Index',axes=F)
y_axis = axis(2,lwd=3)
axis(1,lwd=3,at=seq(0,180,by=20))
p_limits = par("usr");
lines(angle_search$test_angles,abs(angle_search$mean_angle)+20)


highlight_points = c(0,45,90,135);
for (angle in highlight_points) {
    degree_index = which(angle_search$test_angles == angle)
    points(angle_search$test_angles[degree_index],
        angle_search$angle_FAAI[degree_index],col='red',pch=19)
}

points(best_angle, 90-sd(corrected_orientation),col='green',pch=19)

graphics.off()
