rm(list = ls());
source('FA_analysis_lib.R');
library(lattice);
library(geneplotter);
library(Hmisc);

debug = TRUE;

################################################################################
#Result loading
################################################################################
raw_data <- list()

#Rapamycin results
# exp_dirs_rap <- Sys.glob('../../results/rap_src/*/adhesion_props/models/');
# exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)];
# raw_data$rap = load_results(exp_dirs_rap,file.path('intensity.Rdata'));
# 
# exp_dirs_rap <- Sys.glob('../../results/rap_src_control/*/adhesion_props/models/');
# exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)];
# raw_data$rap_ctrl = load_results(exp_dirs_rap,file.path('intensity.Rdata'));
# 
# exp_dirs_rap <- Sys.glob('../../results/rap_src_pax/*/adhesion_props/models/');
# exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)];
# raw_data$rap_pax = load_results(exp_dirs_rap,file.path('intensity.Rdata'));

exp_dirs_rap <- Sys.glob('../../results/rap_fak_no_split/*/adhesion_props/models/');
post_dirs = exp_dirs_rap[regexpr('post', exp_dirs_rap) > -1];
pre_dirs = exp_dirs_rap[regexpr('pre', exp_dirs_rap) > -1];
stopifnot(length(post_dirs) + length(pre_dirs) == length(exp_dirs_rap));

raw_data$rap_fak_pre$intensity = load_results(pre_dirs,file.path('intensity.Rdata'));
raw_data$rap_fak_pre$static_props <- load_data_files(pre_dirs, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

raw_data$rap_fak_post$intensity = load_results(post_dirs,file.path('intensity.Rdata'));
raw_data$rap_fak_post$static_props <- load_data_files(post_dirs, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

print('Done Loading Data')

################################################################################
#Processing
################################################################################

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

        processed$no_filt[[exp_type]] = filter_results(raw_data[[exp_type]]$intensity, 
            min_R_sq = -Inf, max_p_val = Inf, pos_slope=FALSE);
        
        processed$only_signif[[exp_type]] = filter_results(raw_data[[exp_type]]$intensity, 
            min_R_sq = -Inf, max_p_val = 0.05);
    }
    
    dynamic_props[[exp_type]] = gather_general_dynamic_props(raw_data[[exp_type]]$intensity)
    static_props[[exp_type]] = gather_static_props(raw_data[[exp_type]]$static_props) 
}


out_folder = '../../doc/publication/figures'
stop();

################################################################################
#Plotting
################################################################################
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

pre = list()
post = list()

pre$dynamic$death_rate = unlist(lapply(raw_data$rap_fak_pre$intensity, determine_death_rate))
post$dynamic$death_rate = unlist(lapply(raw_data$rap_fak_post$intensity, determine_death_rate))

pre$dynamic$birth_rate = unlist(lapply(raw_data$rap_fak_pre$intensity, determine_birth_rate))
post$dynamic$birth_rate = unlist(lapply(raw_data$rap_fak_post$intensity, determine_birth_rate))

for (i in names(static_props$rap_fak_pre)) {
    if (i == "exp_num") {
        next;
    }
    print(i);
    pre$static[[i]] = tapply(static_props$rap_fak_pre[[i]], static_props$rap_fak_pre$exp_num, mean)
    post$static[[i]] = tapply(static_props$rap_fak_post[[i]], static_props$rap_fak_post$exp_num, mean)
}

for (i in names(dynamic_props$rap_fak_pre)) {
    if (i == "exp_num") {
        next;
    }
    print(i);
    pre$dynamic[[i]] = tapply(dynamic_props$rap_fak_pre[[i]], dynamic_props$rap_fak_pre$exp_num, mean, na.rm=T)
    post$dynamic[[i]] = tapply(dynamic_props$rap_fak_post[[i]], dynamic_props$rap_fak_post$exp_num, mean, na.rm=T)
}

########################################
#RAP-SRC Plotting
########################################

#Lifetime Plots
dir.create(dirname(file.path(out_folder,'rapr_src','rapr_src_rates_nofilt.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'rapr_src','rapr_FAK_lifetimes.svg'))
par(bty='n', mar=c(3,4,0,0))
temp = plot_stage_length_data(rap_src_stage_lengths, type='side_by_side', names=c('Before', 'After'));
graphics.off()

svg(file.path(out_folder,'rapr_src','rapr_pax_lifetimes.svg'))
par(bty='n', mar=c(3,4,0,0))
temp = plot_stage_length_data(rap_pax_stage_lengths, type='side_by_side', names=c('Before', 'After'));
graphics.off()

svg(file.path(out_folder,'rapr_src','rapr_overall_lifetimes.svg'))
plot_max = max(pax_lifetimes$conf_ints, fak_lifetimes$conf_ints) + 4;
par(bty='n', mar=c(2,4,0,0))
barplot(cbind(fak_lifetimes$median_vals, pax_lifetimes$median_vals), 
        beside=TRUE, names=c('FAK', 'Paxillin'), legend=c('Before', 'After'), 
        xlim=c(0,7.5), ylim=c(0,plot_max), ylab='Lifetime (min)')
errbar(c(1.5,2.5,4.5,5.5),c(fak_lifetimes$median_vals, pax_lifetimes$median_vals), 
       c(fak_lifetimes$conf_ints[,2],pax_lifetimes$conf_ints[,2]),
       c(fak_lifetimes$conf_ints[,1],pax_lifetimes$conf_ints[,1]), add=TRUE, cex=1E-20, xlab='', ylab='')

pax_max = max(pax_lifetimes$conf_ints) + 0.75;
upper_left = c(4.5, pax_max+1.5);
lower_right = c(5.5, pax_max);
plot_signif_bracket(upper_left, lower_right, over_text='*')
graphics.off()

#Unfiltered Rate Plotting
pdf(file.path(out_folder,'rapr_src','rapr_FAK_rates_nofilt.pdf'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,1,0))

boxplot_with_points(list(rap_src_pre$assembly$slope, 
                         rap_src_post$assembly$slope),
                    names=c('Before', 'After'), main='Assembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))

boxplot_with_points(list(rap_src_pre$disassembly$slope, 
                         rap_src_post$disassembly$slope),
                    names=c('Before', 'After'), main='Disassembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))

par(bty='n', mar=c(2,4,1.5,0))

boxplot_with_points(list(rap_src_pre_ctrl$assembly$slope, 
                         rap_src_post_ctrl$assembly$slope),
                    names=c('Before', 'After'), main='Assembly Control', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))

boxplot_with_points(list(rap_src_pre_ctrl$disassembly$slope, 
                         rap_src_post_ctrl$disassembly$slope),
                    names=c('Before', 'After'), main='Disassembly Control', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))
graphics.off()

pdf(file.path(out_folder,'rapr_src','rapr_pax_rates_nofilt.pdf'), height=(7/2))
layout(rbind(c(1,2)))
par(bty='n', mar=c(2,4,1,0))

boxplot_with_points(list(rap_pax_pre$assembly$slope, 
                         rap_pax_post$assembly$slope),
                    names=c('Before', 'After'), main='Assembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')), point_cex=0.05)

boxplot_with_points(list(rap_pax_pre$disassembly$slope, 
                         rap_pax_post$disassembly$slope),
                    names=c('Before', 'After'), main='Disassembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')), point_cex=0.05)
graphics.off()


