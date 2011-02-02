source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/presentations/Hahn Lab 8-17-09/figures/R_figures/' 
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

########################################
#Kinetics Figure
#######################################
dir.create(dirname(file.path(out_folder,'kinetics','kinetics_example.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'kinetics','kinetics_example.svg'), width=5, height=5);
layout(rbind(c(1,2),c(3,4)))

exp_one_only = load_results(exp_dirs[[1]],'intensity.Rdata')
exp_one_only = exp_one_only[[1]]

par(bty='n', mar=c(4,4.2,0.5,0))

plot.new()

ad_num = 675
plot_ad_seq(exp_one_only, ad_num, type='overall', 
	phase_lengths=c(exp_one_only$assembly$length[ad_num],exp_one_only$disassembly$length[ad_num]))

plot_ad_seq(exp_one_only, ad_num);
text(4,0.67,pos=3,expression(paste(R^2,' = 0.920')))
text(4,0.57,pos=3,adj=0, 
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',exp_one_only$assembly$slope[ad_num]))))

plot_ad_seq(exp_one_only,ad_num,type='disassembly')
text(4.2,0.33,pos=3, expression(paste(R^2,' = 0.961')))
text(4.2,0.27,pos=3,adj=0,
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',exp_one_only$disassembly$slope[ad_num]))))
graphics.off()


svg(file.path(out_folder,'kinetics','overall_kinetics.svg'), width=5, height=5);
layout(rbind(c(1,2),c(3,3)))
par(bty='n', mar=c(4.1,4.1,1,0.5))
hist(results_nofilt$a$R, main='Assembly', freq=TRUE,
	 xlab='Adjusted R Squared Values', 
	 ylab='# of Focal Adhesions')
hist(results_nofilt$dis$R,main='Disassembly', freq=TRUE,
	 xlab='Adjusted R Squared Values', 
	 ylab='# of Focal Adhesions')

par(bty='n', mar=c(2,4.1,1,0.5))
boxplot_with_points(list(results_filt$a$slope,results_filt$dis$slope), 
		    names=c('Assembly', 'Disassembly'), boxwex=0.6, 
		    ylab=expression(paste('Rate (',min^-1,')',sep='')), point_cex=0.6)
#95% confidence intervals on the mean from Webb 2004
#segments(1.4,0.04,1.4,0.2,lwd=2)
#segments(1.35,0.12,1.45,0.12,lwd=2)
#segments(2.4,0.08,2.4,.088+0.004*2,lwd=2)
#segments(2.35,0.088,2.45,.088,lwd=2)
graphics.off()

########################################
#Spatial Figure
########################################
dir.create(dirname(file.path(out_folder,'spatial','spatial.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'spatial','spatial.svg'), width=7, height=7/2)
par(bty='n',mar=c(4.2,4.1,0.5,0))
layout(rbind(c(1,2)))

breaks_end = ceil(max(c(results_filt$a$edge_dist,results_filt$dis$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);

hist(results_filt$a$edge_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)

hist(results_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)

#par(bty='n',mar=c(4.2,4.1,2,0.4))
#
#max_rate = max(c(results_filt$a$slope, results_filt$d$slope));
#plot(results_filt$a$edge_dist, pch=19, cex=0.4,
#	 results_filt$a$slope,
#	 xlim = c(0,breaks_end),
#         ylim = c(0,max_rate),
#	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
#	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
#
#par(bty='n',mar=c(4.2,4.1,2,0.4))
#
#plot(results_filt$d$edge_dist, pch=19, cex=0.4,
#	 results_filt$d$slope, 
#	 xlim = c(0,breaks_end),
#         ylim = c(0,max_rate),
#	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
#	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
graphics.off()

print('Done with Spacial')

############################################################
#Comparing S178A to Wild-type
############################################################

########################################
#Lifetime Phases
########################################

stage_data <- gather_stage_lengths(corr_results_onlysignif, corr_results_S_onlysignif)
#stage_data_filt <- gather_stage_lengths(results_filt, results_S_filt, debug=TRUE)
bar_lengths = stage_data$bar_lengths
conf_ints = stage_data$conf_ints
conf_ints[1:3,1] = 0.22
conf_ints[4:6,1] = 0.57

err_bars = conf_ints
err_bars[1,3:4] = err_bars[1,3:4]
err_bars[2,3:4] = err_bars[2,3:4] + sum(bar_lengths[1,1])
err_bars[3,3:4] = err_bars[3,3:4] + sum(bar_lengths[1:2,1])
err_bars[4,3:4] = err_bars[4,3:4]
err_bars[5,3:4] = err_bars[5,3:4] + sum(bar_lengths[1,2])
err_bars[6,3:4] = err_bars[6,3:4] + sum(bar_lengths[1:2,2])

dir.create(dirname(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes.svg'), width=6.5, height=5)
par(bty='n',mar=c(2,4,0,0))
barplot(bar_lengths, names=c('Wild-type','S178A'), 
        ylab='Time (min)', width=matrix(0.3,3,2), xlim=c(0,1),
        legend=c('Assembly','Stability','Disassembly'),ylim = c(0,max(err_bars)+1))
#legend('topright',c('Assembly','Stability','Disassembly'), fill=c('gray10','gray','gray'))
errbar(err_bars[,1],err_bars[,2],err_bars[,3],err_bars[,4],add=TRUE,cex=0.0001, xlab='', ylab='')

#Signifcance Bars
bar_length = 1;
sep_from_data = 0.5;

#upside down identification bar
upper_left = c(err_bars[1,1], min(err_bars[1,3],err_bars[4,3]) - sep_from_data);
lower_right = c(err_bars[4,1], min(err_bars[1,3],err_bars[4,3]) - (sep_from_data + bar_length));
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(upper_left[2],lower_right[2],lower_right[2],upper_left[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,lower_right[2]-sep_from_data,"***",cex=1.5)

#right side up identification bar
upper_left = c(err_bars[3,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(err_bars[6,1], max(err_bars[6,4],err_bars[6,4]) + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)

graphics.off()

svg(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes_alt.svg'), width=6.5, height=5)
par(bty='n',mar=c(2.1,4,0,0))
barplot(t(bar_lengths), names=c('Assembly','Stability','Disassembly'), beside=TRUE, xlim=c(0,11),
        ylab='Time (min)', legend=c('Wild-type','S178A'),ylim = c(0,max(conf_ints)+1))

sideways_err_bars = conf_ints;
sideways_err_bars[1,1] = 1.5;
sideways_err_bars[2,1] = 4.5;
sideways_err_bars[3,1] = 7.5;
sideways_err_bars[4,1] = 2.5;
sideways_err_bars[5,1] = 5.5;
sideways_err_bars[6,1] = 8.5;

errbar(sideways_err_bars[,1],sideways_err_bars[,2],sideways_err_bars[,3],sideways_err_bars[,4],add=TRUE,cex=0.0001, xlab='', ylab='');
graphics.off();

print('Done with Lifetime Phase Lengths')

########################################
#Dynamics
########################################
#svg(file.path(out_folder,'S178A','S178A_vs_wild-type.svg'))
dir.create(dirname(file.path(out_folder,'S178A','S178A_vs_wild-type.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'S178A','S178A_vs_wild-type.pdf'),width=5,height=5)
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,0,0))

max_rate = max(c(results_filt$as$slope,results_S_filt$as$slope,
                 results_filt$dis$slope,results_S_filt$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(results_filt$as$slope,results_S_filt$as$slope), 
        names=c('Wild-type','S178A'), 
        colors=c('orange','blue'),
        ylim = c(0,max_rate + 0.012),
        ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
        inc.n.count = FALSE
)

bar_length = .005;
sep_from_data = 0.005;

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"***",cex=1.5)

#Panel Disassembly Rates
boxplot_with_points(list(results_filt$dis$slope,results_S_filt$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate + 0.012),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
    inc.n.count = FALSE
)

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(results_filt$as$edge_dist,results_S_filt$as$edge_dist,
                 results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), na.rm=T)
boxplot_with_points(list(results_filt$as$edge_dist,results_S_filt$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2), 
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')),
    inc.n.count = FALSE
)

bar_length = 1;
sep_from_data = 1;

upper_left = c(1, max_dist + sep_from_data + bar_length);
lower_right = c(2, max_dist + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"**",cex=1.5)

#Panel Death Distances
boxplot_with_points(list(results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2),
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')),
    inc.n.count = FALSE
)

graphics.off()

print('Done with S178A Comparisons')

############################################################
#Distance versus Background corrected correlation
############################################################

all = list()
for (i in 1:length(background_correlation)) {
	all$distances = c(all$distances, background_correlation[[i]]$distances)
	all$correlations = c(all$correlations, background_correlation[[i]]$correlations)
}	
all = as.data.frame(all)

all_S = list()
for (i in 1:length(background_correlation_S)) {
	all_S$distances = c(all_S$distances, background_correlation_S[[i]]$distances)
	all_S$correlations = c(all_S$correlations, background_correlation_S[[i]]$correlations)
}	
all_S = as.data.frame(all_S)

binned_all = bin_corr_data(all, bin_max=200)
binned_all_S = bin_corr_data(all_S, bin_max=200)

dir.create(dirname(file.path(out_folder,'S178A','dist_vs_corr.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'S178A','dist_vs_corr.pdf'))
y_max = max(c(binned_all$upper, binned_all_S$upper), na.rm=TRUE)

plot(binned_all$mids, binned_all$means, ylim=c(0, y_max), 
        xlab=expression(paste('Mean Distance (', mu, 'm)', sep='')), ylab="correlation")
errbar(binned_all$mids, binned_all$means, binned_all$upper, binned_all$lower, add=TRUE, xlab="", ylab="")

#plot(binned_all_S$mids, binned_all_S$means, ylim=c(0, y_max), col='red', add=TRUE)
errbar(binned_all_S$mids, binned_all_S$means, binned_all_S$upper, binned_all_S$lower, add=TRUE,col='red', xlab="", ylab="")
graphics.off()

print('Done with Distance versus Pax Concentration Correlation')
