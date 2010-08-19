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

print('Done Loading Data')

########################################
#Result filtering
########################################

processed = list();
dynamic_props = list();
static_props = list();
for (exp_type in names(raw_data)) {
    for (property in names(raw_data[[exp_type]])) {
        if (property == "static_props") {
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

########################################
#Spacial Figure Preparation
########################################

# In this set of plots we overlay essentially overlay the two plots that are
# made in the above figure into two figures, an assembly and disassembly
# figure. To do this will will first create the histograms and then plots the
# rates versus position. Since the histogram call sets up all the dimensions of
# the plot, we will have to scale the rates to fit into those spaces.

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

########################################
#Kinetics Figure
#######################################
dir.create(dirname(file.path(out_folder,'UCRF_grant','kinetics-spacial.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'UCRF_grant','kinetics-spacial.svg'));
layout(rbind(c(1,2),c(3,3),c(4,5)), heights=c(1,0.75,0.75))

raw_data_wt_one = raw_data$wild_type$intensity[[1]]

par(bty='n', mar=c(4,4.2,1.5,0))

plot.new()
mtext('A',adj=-.19,side=3,line=-0.5,cex=1.5)

ad_num = 1799
ad_num = 675
plot_ad_seq(raw_data_wt_one, ad_num, type='overall', 
	phase_lengths=c(raw_data_wt_one$assembly$length[ad_num],raw_data_wt_one$disassembly$length[ad_num]))
mtext('B',adj=-.19,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(2.1,4.2,1.1,0))
boxplot_with_points(list(processed$only_signif$wild_type$intensity$assembly$slope,
                         processed$only_signif$wild_type$intensity$disassembly$slope),
		    names=c('Assembly', 'Disassembly'), boxwex=0.6, 
		    ylab=expression(paste('Rate (',min^-1,')',sep='')), 
            point_cex=0.6, with.median.props=FALSE)
#95% confidence intervals on the mean from Webb 2004
#segments(1.4,0.04,1.4,0.2,lwd=2)
#segments(1.35,0.12,1.45,0.12,lwd=2)
#segments(2.4,0.08,2.4,.088+0.004*2,lwd=2)
#segments(2.35,0.088,2.45,.088,lwd=2)
mtext('C',adj=-0.085,side=3,line=-0.5,cex=1.5)

########################################
#Spacial Figure
#######################################

par(bty='n',mar=c(4.2,4.1,0.1,5))

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
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)

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
mtext('E',adj=-.225,side=3,line=-1.5,cex=1.5)
graphics.off()

########################################
#Alternate Only Kinetics Figure
#######################################
dir.create(dirname(file.path(out_folder,'UCRF_grant','kinetics-only.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'UCRF_grant','kinetics-only.svg'));
layout(rbind(c(1,2),c(3,3)), heights=c(1,0.75))

raw_data_wt_one = raw_data$wild_type$intensity[[1]]

par(bty='n', mar=c(4,4.2,1.5,0))

plot.new()
mtext('A',adj=-.19,side=3,line=-0.5,cex=1.5)

ad_num = 1799
ad_num = 675
plot_ad_seq(raw_data_wt_one, ad_num, type='overall', 
	phase_lengths=c(raw_data_wt_one$assembly$length[ad_num],raw_data_wt_one$disassembly$length[ad_num]))
mtext('B',adj=-.19,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(2.1,4.2,1.1,0))
boxplot_with_points(list(processed$only_signif$wild_type$intensity$assembly$slope,
                         processed$only_signif$wild_type$intensity$disassembly$slope),
		    names=c('Assembly', 'Disassembly'), boxwex=0.6, 
		    ylab=expression(paste('Rate (',min^-1,')',sep='')), 
            point_cex=0.6, with.median.props=FALSE)
mtext('C',adj=-0.085,side=3,line=-0.5,cex=1.5)
graphics.off();
