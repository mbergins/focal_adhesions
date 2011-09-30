rm(list = ls())
source('FA_analysis_lib.R');
source('bilinear_modeling.R');
library(lattice);
library(geneplotter);
library(Hmisc);

################################################################################
#Result loading
################################################################################
raw_data <- list()
raw_data_base_dir = '../../results/andrei_no_split/adhesion_bleaching/';

exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'WT/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$WT$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal_length5.Rdata'));

exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'dead/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$dead$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal_length5.Rdata'));

exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'sh2/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$sh2$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal_length5.Rdata'));

exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'sh3/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$sh3$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal_length5.Rdata'));

exp_dirs <- Sys.glob(file.path(raw_data_base_dir,'myr/*/adhesion_props/models/'))
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$myr$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal_length5.Rdata'));

print('Done Loading Data')

############################################################
#Processing
############################################################
norm_ad_counts = list()

for (exp_type in names(raw_data)) {
# for (exp_type in "WT") {
    print(paste("Filtering", exp_type));
    
    norm_ad_counts[[exp_type]]$before = matrix(NA,nrow = length(raw_data[[exp_type]]$intensity),
        ncol=1E5);
    norm_ad_counts[[exp_type]]$after = norm_ad_counts[[exp_type]]$before
    for (exp_num in 1:length(raw_data[[exp_type]]$intensity)) {
        this_exp = raw_data[[exp_type]]$intensity[[exp_num]]
        these_props = this_exp$exp_props
        drug_addition_time = these_props$drug_addition_time[1]
        
        ad_present = !is.nan(this_exp$exp_data)
        ad_present_counts = colSums(ad_present);

        mean_ad_present_before = mean(ad_present_counts[1:drug_addition_time]);
        
        ad_present_counts = ad_present_counts/mean_ad_present_before;
        
        ad_present_before = ad_present_counts[1:(drug_addition_time-1)];
        ad_present_after = ad_present_counts[drug_addition_time:length(ad_present_counts)];

        norm_ad_counts[[exp_type]]$before[exp_num,(1E5-drug_addition_time+2):1E5] = ad_present_before
        norm_ad_counts[[exp_type]]$after[exp_num,1:length(ad_present_after)] = ad_present_after
    }
    norm_ad_counts[[exp_type]]$before = 
        norm_ad_counts[[exp_type]]$before[,!is.na(colSums(norm_ad_counts[[exp_type]]$before))]
    norm_ad_counts[[exp_type]]$after = 
        norm_ad_counts[[exp_type]]$after[,!is.na(colSums(norm_ad_counts[[exp_type]]$after))]
}


out_folder = '../../doc/publication/figures/andrei'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

############################################################
# P-value Calculations
############################################################

system('notify-send "done with R"')
stop();

################################################################################
#Plotting
################################################################################

# for (i in 1:11) {
#     if (i == 1) {
#         plot(c(lowess(norm_ad_counts$WT$before[i,])$y,lowess(norm_ad_counts$WT$after[i,])$y),
#             typ='l',ylim=c(0.9,2))
#     } else {
#         lines(c(lowess(norm_ad_counts$WT$before[i,])$y,lowess(norm_ad_counts$WT$after[i,])$y),typ='l')
#     }
# }
# lines(c(lowess(colMeans(norm_ad_counts$WT$before))$y,lowess(colMeans(norm_ad_counts$WT$after))$y),
#     typ='l',col='red',lwd=3)
 
for (exp_type in names(norm_ad_counts)) {
    size = dim(norm_ad_counts[[exp_type]]$before)
    norm_ad_counts[[exp_type]]$before = norm_ad_counts[[exp_type]]$before[,(size[2] - 40):size[2]]
    norm_ad_counts[[exp_type]]$after = norm_ad_counts[[exp_type]]$after[,1:60]
}

svg(file.path(out_folder,'adhesion_counts.svg'),width=4,height=4);
par(bty='n',mar=c(2.6,2.5,0.1,0), mgp=c(1.6,0.5,0),xpd=T)
plot(c(lowess(colMeans(norm_ad_counts$WT$before))$y,lowess(colMeans(norm_ad_counts$WT$after))$y),
    typ='l',col='black',lwd=3,ylim=c(0.8,1.65),xlab='Time',ylab='Normalized Adhesion Count')
lines(c(lowess(colMeans(norm_ad_counts$dead$before))$y,lowess(colMeans(norm_ad_counts$dead$after))$y),
    typ='l',col='blue',lwd=3)
lines(c(lowess(colMeans(norm_ad_counts$myr$before))$y,lowess(colMeans(norm_ad_counts$myr$after))$y),
    typ='l',col='green',lwd=3)
lines(c(lowess(colMeans(norm_ad_counts$sh2$before))$y,lowess(colMeans(norm_ad_counts$sh2$after))$y),
    typ='l',col='purple',lwd=3)
lines(c(lowess(colMeans(norm_ad_counts$sh3$before))$y,lowess(colMeans(norm_ad_counts$sh3$after))$y),
    typ='l',col='red',lwd=3)

legend('topleft',c('WT','Dead','Myr','SH2','SH3'),
    fill=c('black','blue','green','purple','red'), inset=c(0.01,0),bty='n')
graphics.off()

system("notify-send 'done building images'")
