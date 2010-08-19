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
single_props <- list()

exp_dirs <- Sys.glob('../../results/simulation/*/*/models')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

for (dir in exp_dirs) {
    this_data_type = basename(dirname(dirname(dir)))
    raw_data[[this_data_type]] = load_results(dir,file.path('intensity.Rdata'))
}

print('Done Loading Data')

########################################
#Result filtering
########################################
processed = list();

for (exp_type in names(raw_data)) {
    if (debug) {
        print(paste("Filtering", exp_type));
    }
    processed$no_filt[[exp_type]] = filter_results(raw_data[[exp_type]], 
        min_R_sq = -Inf, max_p_val = Inf, pos_slope= FALSE);
    
    processed$only_signif[[exp_type]] = filter_results(raw_data[[exp_type]], 
        min_R_sq = -Inf, max_p_val = 0.05, pos_slope= FALSE);
    
    processed$high_Rsq[[exp_type]] = filter_results(raw_data[[exp_type]], pos_slope= FALSE);
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
#Stationary
########################################

longevity_filt = raw_data$stationary[[1]]$exp_props$longevity >= 25;
ad_sig_filt = raw_data$stationary[[1]]$exp_props$ad_sig[which(longevity_filt)]
mean_area_filt = raw_data$stationary[[1]]$exp_props$mean_area[which(longevity_filt)]

start_y_filt = raw_data$stationary[[1]]$exp_props$start_x[which(longevity_filt)]
y_clust = kmeans(start_y_filt, seq(min(start_y_filt), max(start_y_filt), length=10))

start_x_filt = raw_data$stationary[[1]]$exp_props$start_y[which(longevity_filt)]
x_clust = kmeans(start_x_filt, seq(min(start_x_filt), max(start_x_filt), length=12))

filt_data = data.frame(ad_sig = ad_sig_filt, mean_area = mean_area_filt, x_clusters = x_clust$cluster,
        y_clusters = y_clust$cluster)

dir.create(dirname(file.path(out_folder, 'simulation', 'stationary_hist.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder, 'simulation', 'stationary_results.svg'), height=7*1.1, width=7*1.1);
par(mar=c(4,4,0.4,0), bty='n');
layout(rbind(c(1,2),c(3,4)))

plot.new()
mtext('A',adj=-.2,side=3,line=-1,cex=1.5);

hist_data = hist(raw_data$stationary[[1]]$exp_props$longevity, main='', xlab='Adhesion Longevity')
text(25, 300+hist_data$counts[which(hist_data$mids == 25)], sum(raw_data$stationary[[1]]$exp_props$longevity >= 25, na.rm=TRUE));
mtext('B',adj=-.2,side=3,line=-1,cex=1.5);

par(mar=c(4,4,1.25,0), bty='n');
boxplot(ad_sig ~ x_clusters, data=filt_data, notch=T, ylab='Average Adhesion Intensity', xlab='Column Number',
        names = 4:15)
for (i in 4:15) {
    predicted_int_vals = seq(0.1405, 0.4724,length=15-3);
    print(median(filt_data$ad_sig[filt_data$x_clusters == i - 3]) - predicted_int_vals[i-3])
    segments(i-0.5-3, predicted_int_vals[i-3], i+0.5-3, predicted_int_vals[i-3], col='red')
}
mtext('C',adj=-.2,side=3,line=-0.75,cex=1.5);

boxplot(mean_area ~ y_clusters, data=filt_data, ylab = 'Mean Adhesion Area', xlab='Row Number')
mtext('D',adj=-.2,side=3,line=-0.75,cex=1.5);

graphics.off()

########################################
#Moving
########################################

ad_sizes = c(1,4,5,12,13,24,29,44,49,68);

moving_filtered = list();
moving_names = c("moving_0_5","moving_1", "moving_2", "moving_3", 
                 "moving_4", "moving_5", "moving_6", "moving_7", 
                 "moving_8", "moving_9", "moving_10");
for (i in moving_names) {
    print(i)
    longev_filt = processed$no_filt[[i]]$assembly$longevity >= 25;
    if (any(longev_filt)) {
        this_min_size = min(processed$no_filt[[i]]$assembly$mean_area[longev_filt])

        moving_filtered$min_size = c(moving_filtered$min_size, this_min_size);
        
        moving_filtered$reliability = c(moving_filtered$reliability,
            sum(processed$no_filt[[i]]$assembly$mean_area[longev_filt] <= this_min_size*1.05 &
                processed$no_filt[[i]]$assembly$mean_area[longev_filt] >= this_min_size*0.95))
        
        moving_filtered$number_detected = c(moving_filtered$number_detected,
            length(which(longev_filt)))
    } else {
        moving_filtered$min_size = c(moving_filtered$min_size, 
            NA);

        moving_filtered$number_detected = c(moving_filtered$number_detected,
            0)
    }
}

moving_filtered$radii = c(seq(0.5,5,by=0.5),NA)
moving_filtered$speeds = c(0.5, seq(1,10,by=1))


svg(file.path(out_folder, 'simulation', 'moving_results.svg'), width=7/2, height=7*(3/4));
par(mar=c(4,4,0,0), bty='n');
layout(rbind(1,2), heights=c(0.45,1))

plot.new()
mtext('A',adj=-.3,side=3,line=-1.6,cex=1.5);

plot(moving_filtered$speed, moving_filtered$radii, xlab = 'Movement Speed (pixels/min)', ylab = 'Minimum FA Radii (pixels)')
mtext('B',adj=-.3,side=3,line=-1.25,cex=1.5);
graphics.off()

########################################
#Phases
########################################

svg(file.path(out_folder, 'simulation', 'kinetics_results.svg'));
par(mar=c(4,4,0,0), bty='n');
layout(rbind(c(1,2),c(3,4)))

boxplot(processed$no_filt$phases_10$assembly$slope,
        processed$no_filt$phases_11$assembly$slope, 
        processed$no_filt$phases_12$assembly$slope, 
        processed$no_filt$phases_13$assembly$slope, 
        processed$no_filt$phases_14$assembly$slope, 
        processed$no_filt$phases_15$assembly$slope, 
        processed$no_filt$phases_16$assembly$slope, 
        processed$no_filt$phases_17$assembly$slope, 
        processed$no_filt$phases_18$assembly$slope, 
        processed$no_filt$phases_19$assembly$slope, 
        processed$no_filt$phases_20$assembly$slope,
        ylab = "Assembly Phase Slope", xlab = "Assembly Phase Lengths", names=10:20, notch=T);
for (i in 1:11) {
    phase_name = paste("phases_", i+9, sep='');
    # print(median(processed$no_filt[[phase_name]]$assembly$slope) - log(0.4723663/0.2341950)/(i+8))
    segments(i-0.5,log(0.4723663/0.2341950)/(i+8), i+0.5, log(0.4723663/0.2341950)/(i+8), col='red')
}
mtext('A',adj=-.2,side=3,line=-1.5,cex=1.5);

boxplot(processed$no_filt$phases_10$assembly$length,
        processed$no_filt$phases_11$assembly$length, 
        processed$no_filt$phases_12$assembly$length, 
        processed$no_filt$phases_13$assembly$length, 
        processed$no_filt$phases_14$assembly$length, 
        processed$no_filt$phases_15$assembly$length, 
        processed$no_filt$phases_16$assembly$length, 
        processed$no_filt$phases_17$assembly$length, 
        processed$no_filt$phases_18$assembly$length, 
        processed$no_filt$phases_19$assembly$length, 
        processed$no_filt$phases_20$assembly$length,
        ylab = "Detected Assembly Phase Length", xlab = "Assembly Phase Lengths", names=10:20, notch=T);
for (i in 1:11) {
        segments(i-0.5,i+9, i+0.5, i+9, col='red')
}
mtext('B',adj=-.2,side=3,line=-1.5,cex=1.5);

par(mar=c(4,4,1,0), bty='n');
boxplot(processed$no_filt$phases_10$disassembly$slope,
        processed$no_filt$phases_11$disassembly$slope, 
        processed$no_filt$phases_12$disassembly$slope, 
        processed$no_filt$phases_13$disassembly$slope, 
        processed$no_filt$phases_14$disassembly$slope, 
        processed$no_filt$phases_15$disassembly$slope, 
        processed$no_filt$phases_16$disassembly$slope, 
        processed$no_filt$phases_17$disassembly$slope, 
        processed$no_filt$phases_18$disassembly$slope, 
        processed$no_filt$phases_19$disassembly$slope, 
        processed$no_filt$phases_20$disassembly$slope,
        ylab = "Disassembly Phase Slope", xlab = "Disassembly Phase Lengths", names=10:20, notch=T);
for (i in 1:11) {
        segments(i-0.5,log(0.4723663/0.2341950)/(i+8), i+0.5, log(0.4723663/0.2341950)/(i+8), col='red')
}
mtext('C',adj=-.2,side=3,line=-1,cex=1.5);

boxplot(processed$no_filt$phases_10$disassembly$length,
        processed$no_filt$phases_11$disassembly$length, 
        processed$no_filt$phases_12$disassembly$length, 
        processed$no_filt$phases_13$disassembly$length, 
        processed$no_filt$phases_14$disassembly$length, 
        processed$no_filt$phases_15$disassembly$length, 
        processed$no_filt$phases_16$disassembly$length, 
        processed$no_filt$phases_17$disassembly$length, 
        processed$no_filt$phases_18$disassembly$length, 
        processed$no_filt$phases_19$disassembly$length, 
        processed$no_filt$phases_20$disassembly$length,
        ylab = "Detected Disassembly Phase Length", xlab = "Disassembly Phase Lengths", names=10:20, notch=T);
for (i in 1:11) {
        segments(i-0.5,i+9, i+0.5, i+9, col='red')
}
mtext('D',adj=-.2,side=3,line=-1,cex=1.5);
graphics.off()
