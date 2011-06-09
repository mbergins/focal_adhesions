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
Pax_p = determine_mean_p_value(IA32_Pax$best_FAAI,KD_Pax$best_FAAI)
FAK_p = determine_mean_p_value(IA32_FAK$best_FAAI,KD_FAK$best_FAAI)
Vin_p = determine_mean_p_value(IA32_Vin$best_FAAI,KD_Vin$best_FAAI)

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
par(bty='n',mar=c(2.6,2.6,0.3,0), mgp=c(1.6,0.5,0),xpd=T)

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

###########################################################
# Single Adhesion Stability
###########################################################

source('FA_alignment_search.R')
align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_trial_2/WT_*/*/FA_orientation.Rdata'));
lin_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_1ug_trial_2/WT_*/adhesion_props/single_lin.csv'));
ad_dev_1ug_no_filt = gather_all_single_adhesion_deviances(align_files,lin_files,min.area=-Inf,min.data.points=2);
save(ad_dev_1ug_no_filt,file='../../results/emma/processed_2stdev/Fibro_1ug_trial_2/WT_01/adhesion_props/all_1ug_sing_ad_dev.Rdata')

source('FA_alignment_search.R')
align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_trial_2/WT_*/*/FA_orientation.Rdata'));
lin_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_10ug_trial_2/WT_*/adhesion_props/single_lin.csv'));
ad_dev_10ug_no_filt = gather_all_single_adhesion_deviances(align_files,lin_files,min.area=-Inf,min.data.points=2);
save(ad_dev_10ug_no_filt,file='../../results/emma/processed_2stdev/Fibro_10ug_trial_2/WT_01/adhesion_props/all_10ug_sing_ad_dev.Rdata')
system('notify-send "done with R"')

source('FA_alignment_search.R')
align_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_trial_2/WT_*/*/FA_orientation.Rdata'));
lin_files = Sys.glob(file.path(raw_data_base_dir,'Fibro_100ug_trial_2/WT_*/adhesion_props/single_lin.csv'));
ad_dev_100ug = gather_all_single_adhesion_deviances(align_files,lin_files);
system('notify-send "done with R"')

save(ad_dev_100ug,file='../../results/emma/processed_2stdev/Fibro_100ug_trial_2/WT_01/adhesion_props/all_100ug_single_ad_dev.Rdata');

svg(file.path(out_folder,'FAAI','single_ad_stability_wide.svg'), width=10);
par(bty='n',mar=c(2.8,2.6,.5,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
hist(ad_dev_100ug$mean_dev,xlab='Average Deviance from Start (degrees)',
	ylab='# of Focal Adhesions',main='',axes=F, lwd=3,ylim=c(0,2500))
axis(2,at=seq(0,2500,by=500),lwd=3)
axis(1,lwd=3)
graphics.off()

svg(file.path(out_folder,'FAAI','single_ad_stability_wide_w_count.svg'), width=10);
par(bty='n',mar=c(2.8,2.6,.5,0), mgp=c(1.6,0.5,0),xpd=T,cex=2,lwd=3)
hist(ad_dev_100ug$mean_dev,xlab='Average Deviance from Start (degrees)',
	ylab='# of Focal Adhesions',main='',axes=F, lwd=3,ylim=c(0,2500))
axis(2,at=seq(0,2500,by=500),lwd=3)
axis(1,lwd=3)
text(20,200,paste('n=',dim(ad_dev_100ug)[1]),pos=4)
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

all_spatial$class = floor(all_spatial$dists/20);
concen_summary = list();
for (this_concen in unique(all_spatial$concen)) {
    temp = list();
    sum = 0;
    for (class_num in seq(0,20)) {
        this_class_set = subset(all_spatial,class==class_num & concen == this_concen);
        sum = sum+dim(this_class_set)[1]

        temp$dist = c(temp$dist, mean(this_class_set$dist));
        temp$n_count = c(temp$n_count, length(this_class_set$dist));
        temp$or_mean = c(temp$or_mean,mean(this_class_set$or_diff));
        temp$or_plus = c(temp$or_plus,t.test(this_class_set$or_diff,conf.level=0.95)$conf[2]);
        temp$or_minus = c(temp$or_minus,t.test(this_class_set$or_diff,conf.level=0.95)$conf[1]);
    }
    concen_summary[[this_concen]]=temp;
    print(sum)
}

svg(file.path(out_folder,'spatial','fibro_dist_vs_orientation.svg'),width=8,height=4);
layout(cbind(1,2))
par(bty='n',mar=c(2.7,2.6,0.6,0.2), mgp=c(1.6,0.5,0),xpd=T)
errbar(concen_summary$ug_1$dist, concen_summary$ug_1$or_mean,
    concen_summary$ug_1$or_plus,concen_summary$ug_1$or_minus,
    col='green',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',ylim=c(0,40))
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
    ylab='Adhesion Count')
lines(concen_summary$ug_100$dist, concen_summary$ug_100$n_count)
points(concen_summary$ug_10$dist, concen_summary$ug_10$n_count,col='red')
lines(concen_summary$ug_10$dist, concen_summary$ug_10$n_count,col='red')
points(concen_summary$ug_1$dist, concen_summary$ug_1$n_count,col='green')
lines(concen_summary$ug_1$dist, concen_summary$ug_1$n_count,col='green')
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

tag_spatial$class = floor(tag_spatial$dists/20);
tag_summary = list();
for (this_tag in unique(tag_spatial$tag)) {
    for (this_cell_type in unique(tag_spatial$cell_type)) {
        temp = list();
        sum = 0;
        for (class_num in seq(0,20)) {
            this_class_set = subset(tag_spatial,
                class==class_num & tag == this_tag & cell_type == this_cell_type);
            sum = sum+dim(this_class_set)[1]

            temp$dist = c(temp$dist, mean(this_class_set$dist));
            temp$n_count = c(temp$n_count, length(this_class_set$dist));
            temp$or_mean = c(temp$or_mean,mean(this_class_set$or_diff));
            temp$or_plus = c(temp$or_plus,t.test(this_class_set$or_diff,conf.level=0.95)$conf[2]);
            temp$or_minus = c(temp$or_minus,t.test(this_class_set$or_diff,conf.level=0.95)$conf[1]);
        }
        print(sum)
        tag_summary[[this_tag]][[this_cell_type]]=temp;
    }
}

svg(file.path(out_folder,'spatial','tags_dist_vs_orientation.svg'),width=8,height=12);
layout(rbind(c(1,2),c(3,4),c(5,6)))
for (tag_type in c('Pax','Vin','FAK')) {
    par(bty='n',mar=c(2.7,2.6,3,0.2), mgp=c(1.6,0.5,0),xpd=F)
    errbar(tag_summary[[tag_type]]$KD$dist, tag_summary[[tag_type]]$KD$or_mean,
        tag_summary[[tag_type]]$KD$or_plus,tag_summary[[tag_type]]$KD$or_minus,
        col='red',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',
        ylim=c(0,max(tag_summary[[tag_type]]$KD$or_plus)),xlim=c(0,410))
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

errbar(tag_summary$Vin$KD$dist, tag_summary$Vin$KD$or_mean,
    tag_summary$Vin$KD$or_plus,tag_summary$Vin$KD$or_minus,
    col='red',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',
    ylim=c(0,max(tag_summary$Vin$KD$or_plus)))
lines(lowess(tag_summary$Vin$KD$dist,tag_summary$Vin$KD$or_mean),col='red',lwd=3,lty=3)

errbar(tag_summary$Vin$IA32$dist, tag_summary$Vin$IA32$or_mean,
    tag_summary$Vin$IA32$or_plus,tag_summary$Vin$IA32$or_minus,
    col='red',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',add=T)
lines(lowess(tag_summary$Vin$IA32$dist,tag_summary$Vin$IA32$or_mean),col='red',lwd=3)

errbar(tag_summary$FAK$KD$dist, tag_summary$FAK$KD$or_mean,
    tag_summary$FAK$KD$or_plus,tag_summary$FAK$KD$or_minus,
    col='green',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',
    add=T)
lines(lowess(tag_summary$FAK$KD$dist,tag_summary$FAK$KD$or_mean),col='green',lwd=3,lty=3)

errbar(tag_summary$FAK$IA32$dist, tag_summary$FAK$IA32$or_mean,
    tag_summary$FAK$IA32$or_plus,tag_summary$FAK$IA32$or_minus,
    col='green',xlab='Distance Between Adhesions',ylab='Mean Angle Difference',add=T)
lines(lowess(tag_summary$FAK$IA32$dist,tag_summary$FAK$IA32$or_mean),col='green',lwd=3)

errbar(tag_summary$Pax$KD$dist, tag_summary$Pax$KD$or_mean,
    tag_summary$Pax$KD$or_plus,tag_summary$Pax$KD$or_minus,
    xlab='Distance Between Adhesions',ylab='Mean Angle Difference',
    add=T) 
lines(lowess(tag_summary$Pax$KD$dist,tag_summary$Pax$KD$or_mean),lwd=3,lty=3)

errbar(tag_summary$Pax$IA32$dist, tag_summary$Pax$IA32$or_mean,
    tag_summary$Pax$IA32$or_plus,tag_summary$Pax$IA32$or_minus,
    xlab='Distance Between Adhesions',ylab='Mean Angle Difference',add=T)
lines(lowess(tag_summary$Pax$IA32$dist,tag_summary$Pax$IA32$or_mean),lwd=3)

source('FA_alignment_search.R');
start_time = proc.time()
data_summary = find_dist_overlaps_and_orientations(a,min.overlap=10,output_file=NA);
# end_time = proc.time()
# end_time - start_time

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
