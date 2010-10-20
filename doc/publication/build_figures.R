rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
debug=TRUE;

################################################################################
#Result loading
################################################################################
raw_data <- list()
single_props <- list()

#Wild-type FA
exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

raw_data$wild_type$intensity = load_results(exp_dirs,file.path('intensity.Rdata'));
raw_data$wild_type$static_props <- load_data_files(exp_dirs, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

#S178A Results
exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]

raw_data$S178A$intensity = load_results(exp_dirs_S,file.path('intensity.Rdata'));
raw_data$S178A$static_props <- load_data_files(exp_dirs_S, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

#Wild-type FAK
exp_dirs <- Sys.glob('../../results/FAK/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

raw_data$FAK$intensity = load_results(exp_dirs,file.path('intensity.Rdata'));
raw_data$FAK$static_props <- load_data_files(exp_dirs, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

print('Done Loading Data')

########################################
#Result filtering
########################################

processed = list();
dynamic_props = list();
static_props = list();
for (exp_type in names(raw_data)) {
    for (property in names(raw_data[[exp_type]])) {
        if (property != "intensity") {
            next;
        }
        if (debug) {
            print(paste("Filtering", exp_type, property));
        }

        processed$no_filt[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
            min_R_sq = -Inf, max_p_val = Inf, pos_slope=FALSE);
        
        processed$only_signif[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
            min_R_sq = -Inf, max_p_val = 0.05);
    }
    
    dynamic_props[[exp_type]] = gather_general_dynamic_props(raw_data[[exp_type]]$intensity)
    static_props[[exp_type]] = gather_static_props(raw_data[[exp_type]]$static_props) 
}

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

stop()

################################################################################
#Plotting
################################################################################

filters = lapply(raw_data$wild_type$intensity, produce_rate_filters, min_R_sq=-Inf)

pos_slope = c(); p_val = c(); good_R_sq = c(); all_death = c()
for (set in filters) {
    all_death = c(all_death, set$filter_sets$disassembly$death_status[set$filter_sets$disassembly$low_p_val])
    pos_slope = c(pos_slope, set$filter_sets$disassembly$pos_slope[set$filter_sets$disassembly$low_p_val])
    good_R_sq = c(good_R_sq, set$filter_sets$disassembly$good_R_sq[set$filter_sets$disassembly$low_p_val])
    p_val = c(p_val, set$filter_sets$disassembly$low_p_val)
}
 
all = cbind(pos_slope, good_R_sq, all_death)
counts = vennCounts(all)
vennDiagram(counts)

pos_slope = c(); p_val = c(); good_R_sq = c(); sb = c()
for (set in filters) {
    sb = c(sb, set$filter_sets$assembly$not_split_birth[set$filter_sets$assembly$low_p_val])
    pos_slope = c(pos_slope, set$filter_sets$assembly$pos_slope[set$filter_sets$assembly$low_p_val])
    good_R_sq = c(good_R_sq, set$filter_sets$assembly$good_R_sq[set$filter_sets$assembly$low_p_val])
    p_val = c(p_val, set$filter_sets$assembly$low_p_val)
}
 
all = cbind(pos_slope, good_R_sq, sb)
counts = vennCounts(all)
vennDiagram(counts)

########################################
#Statics Properties
########################################
dir.create(dirname(file.path(out_folder,'statics','statics.svg')), 
    recursive=TRUE, showWarnings=FALSE)
svg(file.path(out_folder,'statics','statics.svg'),height=8)

layout_mat = rbind(c(rep(1,4),rep(2,4),rep(3,4)),
                   c(rep(4,6),rep(5,6)),
                   c(rep(6,6),rep(7,6)),
                   c(rep(0,3), rep(8,6), rep(0, 3)))
layout(layout_mat,heights=c(1,0.65,0.65,0.65))
par(bty='n', mar=c(0,4,1.6,0))

#Place holders for cell images
plot.new()
mtext('A',adj=-.31,side=3,line=0,cex=1.5)

plot.new()
mtext('B',adj=-.31,side=3,line=0,cex=1.5)

plot.new()
mtext('C',adj=-.31,side=3,line=0,cex=1.5)

#Histograms
par(bty='n', mar=c(4,4.2,0,0.1))
area_data = static_props$wild_type$Area[static_props$wild_type$Area < 5];
area_hist = hist(area_data, main="", ylab = "FA Count", 
	 xlab = 'FA Area (\u03BCm\u00B2)');
    
#area_hist_model = lm(log(area_hist$counts) ~area_hist$mids);
#predictions = exp(predict(area_hist_model))
#lines(area_hist$mids,predictions, col='red', lwd=2)

mtext('D',adj=-.2,side=3,line=-0.25,cex=1.5)
hist(static_props$wild_type$ad_sig, main="", ylab = "FA Count", xlab = "Normalized Average Paxillin Intensity")

par(bty='n', mar=c(4,4.2,1.2,0.1))
hist(static_props$wild_type$ax[static_props$wild_type$ax < 8], main="", ylab = "FA Count",  xlab = "Axial Ratio")
hist(static_props$wild_type$cent_dist, main="", ylab = "FA Count",  
	 xlab = expression(paste("Distance from Edge (", symbol("m"), m, ')',sep='')))
hist(dynamic_props$wild_type$longevity, main="", ylab = "FA Count",  xlab = "Longevity (min)")
graphics.off()

svg(file.path(out_folder,'statics','longev_inset.svg'), width=3, height=3/2)
par(bty='n', mar=c(2,2,0.5,0))
hist(dynamic_props$wild_type$longevity[!is.na(dynamic_props$wild_type$longevity) & 
                                       dynamic_props$wild_type$longevity > 20], 
    main="", ylab = "", xlab = "")
graphics.off()

print('Done with Static Properties')

########################################
#Kinetics Figure
#######################################
dir.create(dirname(file.path(out_folder,'kinetics','kinetics.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'kinetics','kinetics.svg'),height=10.5);
layout(rbind(c(1,2),c(3,4),c(5,5)))

raw_data_wt_one = raw_data$wild_type$intensity[[1]]

par(bty='n', mar=c(4,4.2,1.5,0))

plot.new()
mtext('A',adj=-.19,side=3,line=-0.5,cex=1.5)

ad_num = 1799
ad_num = 675
plot_ad_seq(raw_data_wt_one, ad_num, type='overall', 
	phase_lengths=c(raw_data_wt_one$assembly$length[ad_num],raw_data_wt_one$disassembly$length[ad_num]))
mtext('B',adj=-.19,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(4,4.2,4,0))
plot_ad_seq(raw_data_wt_one, ad_num, main = 'Assembly');
limits = par("usr");
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3,
     substitute(paste('R' ^2,' = ', x,sep=''), list(x=sprintf('%.03f',raw_data_wt_one$assembly$R_sq[ad_num]))))
text(3,(limits[4]-limits[3])*0.75+limits[3],pos=3,
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',raw_data_wt_one$assembly$slope[ad_num]))))
mtext('C',adj=-.19,side=3,line=-0.5,cex=1.5)

plot_ad_seq(raw_data_wt_one,ad_num,type='disassembly', main = 'Disassembly')
limits = par("usr");
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3,
     substitute(paste('R' ^2,' = ', x,sep=''), list(x=sprintf('%.03f',raw_data_wt_one$disassembly$R_sq[ad_num]))))
text(3,(limits[4]-limits[3])*0.75+limits[3],pos=3,
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',raw_data_wt_one$disassembly$slope[ad_num]))))
mtext('D',adj=-.19,side=3,line=-0.5,cex=1.5)
par(bty='n', mar=c(2.1,4.2,1.1,0))

boxplot_with_points(list(processed$only_signif$wild_type$intensity$assembly$slope,
                         processed$only_signif$wild_type$intensity$disassembly$slope),
		    names=c('Assembly', 'Disassembly'), boxwex=0.6, 
		    ylab=expression(paste('Rate (',min^-1,')',sep='')), 
            point_cex=0.6, with.p.vals=FALSE)
#95% confidence intervals on the mean from Webb 2004
#segments(1.4,0.04,1.4,0.2,lwd=2)
#segments(1.35,0.12,1.45,0.12,lwd=2)
#segments(2.4,0.08,2.4,.088+0.004*2,lwd=2)
#segments(2.35,0.088,2.45,.088,lwd=2)
mtext('E',adj=-0.085,side=3,line=-0.5,cex=1.5)
graphics.off()


####################
#Confidence Intervals
####################

library(boot)
assembly_conf = boot(processed$only_signif$wild_type$intensity$assembly$slope, 
    function(values,indexes) median(values[indexes],na.rm=T), 10000);
disassembly_conf = boot(processed$only_signif$wild_type$intensity$disassembly$slope, 
    function(values,indexes) median(values[indexes],na.rm=T), 10000);

mean(processed$only_signif$wild_type$intensity$assembly$slope, na.rm=T);
sd(processed$only_signif$wild_type$intensity$assembly$slope, na.rm=T);
mean(processed$only_signif$wild_type$intensity$disassembly$slope, na.rm=T);
sd(processed$only_signif$wild_type$intensity$disassembly$slope, na.rm=T);

####################
#Supplemental
####################
dir.create(dirname(file.path(out_folder,'supplemental','R_squared.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'supplemental','R_squared.svg'), width=7, height=7)
layout(rbind(c(1,2),c(3,4)))

wt_R_sq_data = processed$no_filt$wild_type$intensity;
S178A_R_sq_data = processed$no_filt$S178A$intensity;

par(bty='n', mar=c(4,4.2,2,0))

hist(wt_R_sq_data$assembly$R_sq, main='Wild-type Assembly', freq=TRUE,
	 xlab=substitute(paste('Adjusted R' ^2,' Values (n=', x, ')', sep=''), 
                     list(x=length(wt_R_sq_data$assembly$R_sq))), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(wt_R_sq_data$assem$R)
segments(sorted_r_vals[floor(length(sorted_r_vals)/2)],0, 
         sorted_r_vals[floor(length(sorted_r_vals)/2)],plot_dims[4], col='red', lwd=2)
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

hist(wt_R_sq_data$dis$R,main='Wild-type Disassembly', freq=TRUE,
	 xlab=substitute(paste('Adjusted R' ^2,' Values (n=', x, ')', sep=''), 
                     list(x=length(wt_R_sq_data$disassembly$R_sq))), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(wt_R_sq_data$dis$R)
segments(sorted_r_vals[floor(length(sorted_r_vals)/2)],0, 
         sorted_r_vals[floor(length(sorted_r_vals)/2)],plot_dims[4], col='red', lwd=2)
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(4,4.2,2,0))

hist(S178A_R_sq_data$assembly$R_sq, main='S178A Assembly', freq=TRUE,
	 xlab=substitute(paste('Adjusted R' ^2,' Values (n=', x, ')', sep=''), 
                     list(x=length(S178A_R_sq_data$assembly$R_sq))), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(S178A_R_sq_data$assem$R)
segments(sorted_r_vals[floor(length(sorted_r_vals)/2)],0, 
         sorted_r_vals[floor(length(sorted_r_vals)/2)],plot_dims[4], col='red', lwd=2)
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5)

hist(S178A_R_sq_data$dis$R,main='S178A Disassembly', freq=TRUE,
	 xlab=expression(paste('Adjusted R' ^2,' Values (n=2046)', sep='')), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(S178A_R_sq_data$dis$R)
segments(sorted_r_vals[floor(length(sorted_r_vals)/2)],0, 
         sorted_r_vals[floor(length(sorted_r_vals)/2)],plot_dims[4], col='red', lwd=2)
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)

graphics.off()

print('Done with Kinetics')

########################################
#Spacial Figure
########################################

# In this set of plots we overlay the two plots that are made in the above
# figure into two figures, an assembly and disassembly figure. To do this will
# will first create the histograms and then plots the rates versus position.
# Since the histogram call sets up all the dimensions of the plot, we will have
# to scale the rates to fit into those spaces.

# First we need to collect a few properties of the non-overlayed plots,
# specifically, how many and where the tick marks are placed on the rate
# plotting
wt_only_signif = processed$only_signif$wild_type$intensity;

# There are a few adhesions for which we don't know the location of the closest
# cell edge during birth or death due to the cell edge being outside the
# microscope view
wt_na_filt_assembly = subset(wt_only_signif$assembly, ! is.na(edge_dist))
wt_na_filt_disassembly = subset(wt_only_signif$disassembly, ! is.na(edge_dist))

breaks_end = ceil(max(c(wt_na_filt_assembly$edge_dist,wt_na_filt_disassembly$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);
max_rate = max(c(wt_na_filt_assembly$slope, wt_na_filt_disassembly$slope));

#Draw each of the plots to find out the appropriate axis tick positions
plot(wt_na_filt_assembly$edge_dist, pch=19, cex=0.5,
	 wt_na_filt_assembly$slope,
	 xlim = c(0,breaks_end),
     ylim = c(0,max_rate))
assembly_axis_ticks = axTicks(2);
graphics.off()

plot(wt_na_filt_disassembly$edge_dist, pch=19, cex = 0.5,
	 wt_na_filt_disassembly$slope, 
	 xlim = c(0,breaks_end),
     ylim = c(0,max_rate))
disassembly_axis_ticks = axTicks(2);
graphics.off()


# With the preliminary data collected, we continue into the primary plotting
# function
dir.create(dirname(file.path(out_folder,'spatial','spatial.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'spatial','spatial.svg'),height=4, width=10)
par(bty='n',mar=c(4.2,4.1,0.1,5))
layout(rbind(c(1,2)))

##################
#Assembly Plot
##################
hist(wt_na_filt_assembly$edge_dist, 
     xlab=paste('Distance from Edge at Birth (\u03BCm) n=',length(wt_na_filt_assembly$edge_dist), sep=''), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
assembly_plot_dims = par('usr')

# we adjust the scaling factor slightly down, so that all the points are
# included in the plot size, otherwise the highest point gets cutoff
scaled_assembly_rates = wt_na_filt_assembly$slope*((assembly_plot_dims[4]*0.97)/max_rate)

points(wt_na_filt_assembly$edge_dist, scaled_assembly_rates,
       pch=19, cex=0.05, col='darkgreen',
	   ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))

# Now we have to place the axis ticks in the proper location, so we also scale
# the axis ticks, using the same number of marks as would be used in the
# standard plots
axis(4, at = assembly_axis_ticks*((assembly_plot_dims[4]*0.97)/max_rate), 
     labels=assembly_axis_ticks, col='darkgreen', col.axis='darkgreen')

mtext(expression(paste('Assembly Rate (',min^-1,')',sep='')),side=4,line=3, col='darkgreen');
mtext('A',adj=-.22,side=3,line=-1.5,cex=1.5)

##################
#Disassembly Plot
##################
par(bty='n',mar=c(4.2,4.1,0.1,4))
hist(wt_na_filt_disassembly$edge_dist, 
     xlab=paste('Distance from Edge at Death (\u03BCm) n=',length(wt_na_filt_disassembly$edge_dist), sep=''), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
disassembly_plot_dims = par('usr')

# we adjust the scaling factor slightly down, so that all the points are
# included in the plot size, otherwise the highest point gets cutoff
scaled_disassembly_rates = wt_na_filt_disassembly$slope*((disassembly_plot_dims[4]*0.97)/max_rate)

points(wt_na_filt_disassembly$edge_dist, scaled_disassembly_rates,
       pch=19, cex=0.05, col='red',
	   ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))

# Now we have to place the axis ticks in the proper location, so we also scale
# the axis ticks, using the same number of marks as would be used in the
# standard plots
axis(4, at = disassembly_axis_ticks*((disassembly_plot_dims[4]*0.97)/max_rate), 
     labels=disassembly_axis_ticks, col='red', col.axis='red')

mtext(expression(paste('Disassembly Rate (',min^-1,')',sep='')),side=4,line=3, col='red');
mtext('B',adj=-.22,side=3,line=-1.5,cex=1.5)
graphics.off()

####################
#Supplemental
####################
svg(file.path(out_folder,'supplemental','birth_vs_death_pos.svg'))
par(bty='n',mar=c(4.2,4.1,2,0.2))

plot(wt_no_filt$j$birth_dist, wt_no_filt$j$death_dist, 
     xlab='Distance from Edge at Birth (\u03BCm)',
     ylab='Distance from Edge at Death (\u03BCm)', 
     pch=20, cex=0.75)

birth_vs_death_model <- lm(death_dist ~ birth_dist, data=wt_no_filt$joint)
abline(birth_vs_death_model, col='green', lwd = 3)
model_summary <- summary(birth_vs_death_model)
x_data <- data.frame(birth_dist = seq(min(wt_no_filt$j$b, na.rm=T),max(wt_no_filt$j$b, na.rm=T),by=0.01))
line_conf = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, line_conf[,2], col='red', lty=2, lwd = 3)
lines(x_data$birth_dist, line_conf[,3], col='red', lty=2, lwd = 3)

segments(0,0,max(wt_no_filt$j$b, na.rm=T), max(wt_no_filt$j$b, na.rm=T), col='blue', lty=4, lwd = 3)

graphics.off()

####################
# Birth/Death position versus R squared
####################
range_text = c();

assembly_ranges = rep(NA, length(wt_only_signif$assembly$edge_dist));
disassembly_ranges = rep(NA, length(wt_only_signif$disassembly$edge_dist));
for (i in 2:length(these_breaks)) {
    this_filter = wt_only_signif$assembly$edge_dist < these_breaks[i] & 
                  wt_only_signif$assembly$edge_dist >= these_breaks[i - 1]

    assembly_ranges[this_filter] = mean(c(these_breaks[i], these_breaks[i-1]));
    
    this_filter = wt_only_signif$disassembly$edge_dist < these_breaks[i] & 
                  wt_only_signif$disassembly$edge_dist >= these_breaks[i - 1]

    disassembly_ranges[this_filter] = mean(c(these_breaks[i], these_breaks[i-1]));

    range_text = c(range_text, paste(these_breaks[i-1], "-", these_breaks[i]));
}

svg(file.path(out_folder,'supplemental','pos_vs_R_sq.svg'))
layout(c(1,2))

par(bty='n',mar=c(4.2,4.1,2,0.2))
boxplot(wt_only_signif$assembly$R_sq ~ assembly_ranges, names = range_text, las=2, ylab='Adjusted R\u00B2 Values')
mtext("Position at Birth Range (\u03BCm)", side=1, line=4);

par(bty='n',mar=c(5.2,4.1,2,0.2))
boxplot(wt_only_signif$disassembly$R_sq ~ disassembly_ranges, names = range_text, las=2, ylab='Adjusted R\u00B2 Values')
mtext("Position at Death Range (\u03BCm)", side=1, line=4);

graphics.off();
print('Done with Spatial')

############################################################
#Comparing S178A to Wild-type
############################################################

########################################
# Statics Comparisons
########################################
svg(file.path(out_folder,'S178A','statics_comparisons.svg'), height=3.5)
layout(cbind(1,2))
par(bty='n', mar=c(2,4,0,0))

max_area = max(c(static_props$wild_type$Area, static_props$S178A$Area));

boxplot_with_points(list(static_props$wild_type$Area, static_props$S178A$Area), 
    names=c('Wild-type','S178A'), ylab='FA Area (\u03BCm\u00B2)', range=0, inc.points=FALSE,
    inc.n.counts=FALSE)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

boxplot_with_points(list(static_props$wild_type$ax, static_props$S178A$ax), 
    names=c('Wild-type','S178A'), ylab='Axial Ratio', range=0, inc.points=FALSE,
    inc.n.counts=FALSE)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)	    
graphics.off()

########################################
#Lifetime Phases
########################################
stage_data <- gather_stage_lengths(processed$only_signif$wild_type$intensity, 
                                   processed$only_signif$S178A$intensity)

dir.create(dirname(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes.svg'))
par(bty='n',mar=c(2,4,0,0))

err_bars = plot_stage_length_data(stage_data, top_gap=3.5, names=c('Wild-type','S178A'))

#Signifcance Bars
bar_length = 1;
sep_from_data = 0.5;

#upside down identification bar
upper_left = c(err_bars[1,1], min(err_bars[1,3],err_bars[4,3]) - sep_from_data);
lower_right = c(err_bars[4,1], min(err_bars[1,3],err_bars[4,3]) - (sep_from_data + bar_length));
plot_signif_bracket(upper_left,lower_right, orientation='upside_down', over_text = "**", text_x_adj=-0.005)

#right side up identification bar
upper_left = c(err_bars[3,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(err_bars[6,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data);
plot_signif_bracket(upper_left,lower_right, over_text = "*", text_x_adj=-0.005)

graphics.off()

svg(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes_alt.svg'))
par(bty='n',mar=c(2.1,4,0,0))

sideways_err_bars = plot_stage_length_data(stage_data, type='side_by_side', top_gap=1.5, 
    names=c('Wild-type','S178A'), sideways_high_xlim=11);

#Signifcance Bars
bar_length = 0.4;
sep_from_data = 0.2;

#assembly significance bar
upper_left = c(sideways_err_bars[1,1], max(sideways_err_bars[1,4],sideways_err_bars[4,4]) + sep_from_data + bar_length);
lower_right = c(sideways_err_bars[4,1], max(sideways_err_bars[1,4],sideways_err_bars[4,4]) + sep_from_data);
plot_signif_bracket(upper_left,lower_right, over_text = "**")

#disassembly significance bar
upper_left = c(sideways_err_bars[3,1], max(sideways_err_bars[3,4],sideways_err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(sideways_err_bars[6,1], max(sideways_err_bars[3,4],sideways_err_bars[6,4]) + sep_from_data);
plot_signif_bracket(upper_left,lower_right, over_text = "*")

graphics.off();

print('Done with Lifetime Phase Lengths')

########################################
#Dynamics
########################################
dir.create(dirname(file.path(out_folder,'S178A','S178A_vs_wild-type.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'S178A','S178A_vs_wild-type.svg'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4.2,0,0))

wt_only_signif = processed$only_signif$wild_type$intensity;
S178A_only_signif = processed$only_signif$S178A$intensity;

max_rate = max(c(wt_only_signif$as$slope,S178A_only_signif$as$slope,
                 wt_only_signif$dis$slope,S178A_only_signif$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(wt_only_signif$as$slope,S178A_only_signif$as$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate),
    ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
    inc.points=FALSE,
    p.vals.pos = c(0.5,0.9)
)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

#Panel Disassembly Rates
boxplot_with_points(list(wt_only_signif$dis$slope,S178A_only_signif$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
    inc.points = FALSE,
    p.vals.pos = c(0.5,0.9)
)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(wt_only_signif$as$edge_dist,S178A_only_signif$as$edge_dist,
                 wt_only_signif$dis$edge_dist,S178A_only_signif$dis$edge_dist), na.rm=T)
boxplot_with_points(list(wt_only_signif$as$edge_dist,S178A_only_signif$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist), 
    colors=c('orange','blue'),
    ylab='Distance from Edge at Birth (\u03BCm)',
    inc.points = FALSE,
    p.vals.pos = c(0.5,0.9)
)
mtext('C',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Death Distances
boxplot_with_points(list(wt_only_signif$dis$edge_dist,S178A_only_signif$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist),
    colors=c('orange','blue'),
    ylab='Distance from Edge at Death (\u03BCm)', 
    inc.points = FALSE,
    p.vals.pos = c(0.5,0.9)
)
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)	    

graphics.off()

interval_str = list();
interval_str$assembly = print_ratio_conf_string(wt_only_signif$as$slope,S178A_only_signif$as$slope);
interval_str$disassembly = print_ratio_conf_string(wt_only_signif$dis$slope,S178A_only_signif$dis$slope);
interval_str$edge_birth = print_ratio_conf_string(wt_only_signif$as$edge_dist,S178A_only_signif$as$edge_dist);
interval_str$edge_death = print_ratio_conf_string(wt_only_signif$dis$edge_dist,S178A_only_signif$dis$edge_dist);

########################################
#Dynamics - Alternative with significance bar
########################################
dir.create(dirname(file.path(out_folder,'S178A','S178A_vs_wild-type.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'S178A','S178A_vs_wild-type_sig_bar.svg'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,0,0))

wt_only_signif = processed$only_signif$wild_type$intensity;
S178A_only_signif = processed$only_signif$S178A$intensity;

max_rate = max(c(wt_only_signif$as$slope,S178A_only_signif$as$slope,
                 wt_only_signif$dis$slope,S178A_only_signif$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(wt_only_signif$as$slope,S178A_only_signif$as$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    # we increase the the max_rate variable to make space for the significance
    # bar
    ylim = c(0,max_rate*1.1),
    ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
    inc.points=FALSE,
    with.p.vals=FALSE
)

bar_length = max_rate*0.05;
sep_from_data = max_rate*0.025;

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
plot_signif_bracket(upper_left, lower_right,over_text='*');
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

#Panel Disassembly Rates
boxplot_with_points(list(wt_only_signif$dis$slope,S178A_only_signif$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    # we increase the the max_rate variable to make space for the significance
    # bar
    ylim = c(0,max_rate*1.1),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
    inc.points = FALSE,
    with.p.vals=FALSE
)

bar_length = max_rate*0.05;
sep_from_data = max_rate*0.025;

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
plot_signif_bracket(upper_left, lower_right,over_text='*');
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(wt_only_signif$as$edge_dist,S178A_only_signif$as$edge_dist,
                 wt_only_signif$dis$edge_dist,S178A_only_signif$dis$edge_dist), na.rm=T)
boxplot_with_points(list(wt_only_signif$as$edge_dist,S178A_only_signif$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    # we increase the the max_rate variable to make space for the significance
    # bar
    ylim = c(0,max_dist*1.1),
    colors=c('orange','blue'),
    ylab='Distance from Edge at Birth (\u03BCm)',
    inc.points = FALSE,
    with.p.vals=FALSE
)

bar_length = max_dist*0.05;
sep_from_data = max_dist*0.025;

upper_left = c(1, max_dist + sep_from_data + bar_length);
lower_right = c(2, max_dist + sep_from_data);
plot_signif_bracket(upper_left, lower_right,over_text='*');
mtext('C',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Death Distances
boxplot_with_points(list(wt_only_signif$dis$edge_dist,S178A_only_signif$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    # we increase the the max_rate variable to make space for the significance
    # bar
    ylim = c(0,max_dist*1.1),
    colors=c('orange','blue'),
    ylab='Distance from Edge at Death (\u03BCm)', 
    inc.points = FALSE,
    with.p.vals=FALSE
)
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)	    

graphics.off()
print('Done with S178A Comparisons')

#############################################################
#FAK 
#############################################################

########################################
#Dynamics
########################################

# dir.create(file.path(out_folder,'FAK'), recursive=TRUE, showWarnings=FALSE);
# svg(file.path(out_folder,'FAK','FAK_vs_wild-type.svg'))
# layout(rbind(c(1,2),c(3,4)))
# par(bty='n', mar=c(2,4.2,1,0))
# 
# boxplot(list(static_props$wild_type$Area, static_props$FAK$Area), 
#     names=c('Wild-type','FAK'), ylab='FA Area (\u03BCm\u00B2)', range=0)
# mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    
# 
# (mean(static_props$wild_type$Area) - mean(static_props$FAK$Area))/mean(static_props$wild_type$Area)
# 
# boxplot(list(static_props$wild_type$ax, static_props$FAK$ax), 
#     names=c('Wild-type','FAK'), ylab='Axial Ratio', range = 0 )
# mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)	    
# 
# (mean(static_props$wild_type$ax) - mean(static_props$FAK$ax))/mean(static_props$wild_type$ax)

dir.create(file.path(out_folder,'FAK'), recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'FAK','FAK_vs_wild-type_dynamics.svg'), height=3.5, width=8)
layout(cbind(1,2))
par(bty='n', mar=c(2,4.2,0.2,0))

wt_only_signif = processed$only_signif$wild_type$intensity;
FAK_only_signif = processed$only_signif$FAK$intensity;

max_rate = max(c(wt_only_signif$as$slope,FAK_only_signif$as$slope,
                 wt_only_signif$dis$slope,FAK_only_signif$dis$slope));

#Panel Assembly Rates
par(bty='n', mar=c(2,4.2,0.2,-1))
boxplot_with_points(list(wt_only_signif$as$slope,FAK_only_signif$as$slope), 
    names=c('Paxillin','FAK'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate),
    ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
    inc.points=FALSE,
    p.vals.pos = c(0.5,0.9)
)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

#Panel Disassembly Rates
boxplot_with_points(list(wt_only_signif$dis$slope,FAK_only_signif$dis$slope), 
    names=c('Paxillin','FAK'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
    inc.points = FALSE,
    p.vals.pos = c(0.5,0.9)
)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

graphics.off()
print('Done with FAK Comparisons')
