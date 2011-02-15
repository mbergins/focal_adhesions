################################################################################
# FA_analysis_lib.R: various functions used to find and plot the linear regions
#  and associated data from the focal adhesion identification/analysis programs
################################################################################

########################################
#Plotting Functions
########################################
plot_ad_seq <- function (results,index,type='assembly',log.trans = TRUE, time.spacing=1, ...) {
	ad_seq = as.vector(results$exp_data[index,])
	ad_seq = t(ad_seq[!(is.nan(ad_seq))])
	
	if (type == 'assembly') {
		this_ad_seq = ad_seq[1:results$assembly$length[index]];
		stopifnot(results$assembly$length[index] == length(this_ad_seq))
		if (log.trans) {
			this_ad_seq = log(this_ad_seq/this_ad_seq[1]);
		} else {
			this_ad_seq = this_ad_seq/this_ad_seq[1];
		}
		
		x = c(1,results$assembly$length[index]);
		y = c(results$assembly$slope[index]*x[1] + results$assembly$inter[index],
			  results$assembly$slope[index]*x[2] + results$assembly$inter[index])
		
		plot(x[1]:x[2],this_ad_seq,xlab='Time (minutes)',ylab='ln(Intensity/Initial Intensity)',
				 ylim=c(min(this_ad_seq,y),max(this_ad_seq,y)), ...)
		
		lines(x,y,col='darkgreen',lwd = 3)
	}

	if (type == 'disassembly') {
		this_ad_seq = ad_seq[(length(ad_seq) - results$disassembly$length[index] + 1) : length(ad_seq)];
		stopifnot(results$disassembly$length[index] == length(this_ad_seq))
		if (log.trans) {
			this_ad_seq = log(this_ad_seq[1]/this_ad_seq);
		} else {
			this_ad_seq = this_ad_seq[1]/this_ad_seq;
		}
		
		x = c(length(ad_seq) - results$disassembly$length[index] - 1,length(ad_seq));
		y = c(results$disassembly$slope[index]*x[1] + results$disassembly$inter[index],
		   	  results$disassembly$slope[index]*x[2] + results$disassembly$inter[index])
		
		x = c(1,results$disassembly$length[index])
		
		plot(x[1]:x[2],
			 this_ad_seq, xlab='Time (minutes)', ylab='ln(Initial Intensity/Intensity)',
			 ylim=c(min(this_ad_seq,y),max(this_ad_seq,y)), ...)
		
		lines(x,y,col='red',lwd = 3)
	}
	
	if (type == 'overall') {
		x = c(0,results$assembly$length[index]);
		y = c(results$assembly$slope[index]*x[1] + results$assembly$inter[index],
			  results$assembly$slope[index]*x[2] + results$assembly$inter[index]);
		
		plot(0:(length(ad_seq)-1), ad_seq, xlab='Time (minutes)', ylab='Normalized Intensity',type="o");
		
        phase_lengths = c(results$assembly$length[index], results$disassembly$length[index]);

        if (all(! is.na(phase_lengths))) {
            lowess_points = lowess(0:(length(ad_seq)-1), ad_seq, f=1/3)

            lines(lowess_points$x[0:(phase_lengths[1])], lowess_points$y[0:(phase_lengths[1])], col='darkgreen',lwd = 3)
            lines(lowess_points$x[phase_lengths[1]:(length(lowess_points$x) - phase_lengths[2] + 1)], 
                  lowess_points$y[phase_lengths[1]:(length(lowess_points$x) - phase_lengths[2] + 1)],  
                  col='yellow',lwd = 3)
            lines(lowess_points$x[(length(lowess_points$x) - phase_lengths[2] + 1):(length(lowess_points$x))], 
                  lowess_points$y[(length(lowess_points$x) - phase_lengths[2] + 1):(length(lowess_points$x))], 
                  col='red',lwd = 3)
        } else {	
            lines(lowess(0:(length(ad_seq)-1), ad_seq, f=1/3), col='black',lwd = 3)
        }
	}
}

plot_overall_residuals <- function(results,dir,file='overall_residual_plot.pdf',window = 0.1) {
	resid = gather_exp_residuals(results)
	
	resid_win = gather_exp_win_residuals(resid,window=window)
	
	library(Hmisc)
	pdf(file = file.path(dir,file),width=12)
	par(mfcol=c(1,2))
	errbar(resid_win$x$assembly, resid_win$y$assembly, 
		   resid_win$y$assembly - resid_win$err$assembly, resid_win$y$assembly + resid_win$err$assembly, 
		   xlab='Scaled Linear Fit Position',ylab='Residuals')
	lines(c(0,1),c(0,0))
	
	errbar(resid_win$x$disassembly, resid_win$y$disassembly, 
		   resid_win$y$disassembly - resid_win$err$disassembly, resid_win$y$disassembly + resid_win$err$disassembly, 
		   xlab='Scaled Linear Fit Position',ylab='Residuals')
	lines(c(0,1),c(0,0))
	dev.off()
}

gather_exp_residuals <- function(results, min_R_sq = NA) {
	resid = list()
	for (i in 1:length(results)) {
		res = results[[i]]
		for (j in 1:length(res$assembly$R_sq)) {
			if (is.na(res$assembly$R_sq[j])) {	
				next
			}
			if (is.numeric(min_R_sq) & (res$assembly$R_sq[j] < min_R_sq)) {
				next
			}

			resid_list = res$assembly$residual[j][[1]]
			x = list(seq(0,1,1/(length(resid_list)-1)))
			resid$y$assembly = c(resid$y$assembly,resid_list)
			resid$x$assembly = c(resid$x$assembly,x)
		}
		for (j in 1:length(res$disassembly$R_sq)) {
			if (is.na(res$disassembly$R_sq[j])) {
				next
			}
			if (is.numeric(min_R_sq) & (res$disassembly$R_sq[j] < min_R_sq)) {
				next
			}

			resid_list = res$disassembly$residual[j][[1]]
			x = list(seq(0,1,1/(length(resid_list)-1)))
			resid$y$disassembly = c(resid$y$disassembly,resid_list)
			resid$x$disassembly = c(resid$x$disassembly,x)
		}
	}
	resid
}

gather_exp_win_residuals <- function(resid, window) {
	resid_win = list()
	for (i in seq(window,1,by=window)) {
		temp = c()
		for (j in 1:length(resid$x$assembly)) {
			this_resid_x = resid$x$assembly[[j]]
			for (k in 1:length(this_resid_x)) {
				if (this_resid_x[[k]] <= i & this_resid_x[[k]] >= i - window) {
					temp = c(temp,resid$y$assembly[[j]])
				}
			}
		}
		resid_win$x$assembly   = c(resid_win$x$assembly,i-(window/2))
		resid_win$y$assembly   = c(resid_win$y$assembly,mean(temp))
		resid_win$err$assembly = c(resid_win$err$assembly,sd(temp))
		
		temp = c()
		for (j in 1:length(resid$x$disassembly)) {
			this_resid_x = resid$x$disassembly[[j]]
			for (k in 1:length(this_resid_x)) {
				if (this_resid_x[[k]] <= i & this_resid_x[[k]] >= i - window) {
					temp = c(temp,resid$y$disassembly[[j]])
				}
			}
		}
		resid_win$x$disassembly   = c(resid_win$x$disassembly,i-(window/2))
		resid_win$y$disassembly   = c(resid_win$y$disassembly,mean(temp))
		resid_win$err$disassembly = c(resid_win$err$disassembly,sd(temp))
	}
	resid_win
}

boxplot_with_points <- function(data,names=seq(1,length(data),by=1),
    colors=c('darkgreen','red','yellow','blue','pink','cyan','gray','orange','brown','purple'),
    notch=F,  range=1.5, inc.n.counts = TRUE, inc.points = TRUE, pch=20,
    na.omit = TRUE, point.cex=0.2, return_output = FALSE, with.p.values = TRUE, bootstrap.rep = 10000,
    p.values.pos = c(0.5,0.9), p.values.color = 'blue', ...) {
	
    if (any(is.null(data))) {
        print(paste("The data in position", which(is.null(data)), "is null"))
        stop()
    }

    if (na.omit) {
        data = lapply(data,na.omit)
    }
	if (inc.n.counts) {
		for (i in 1:length(data)) {
            names[i] = paste(names[i], ' (n=', length(data[[i]]), ')', sep ='');
		}
	}
	
	box.data = boxplot(data,notch = notch,names = names,varwidth=T,range = range,pch=19,cex=0.25,...)

    plot_dims = par("usr");
	if (inc.points) {
		for (i in 1:length(data)) {
			this_data = data[[i]]
			temp_data = this_data[this_data >= box.data$stat[1,i] & this_data <= box.data$stat[5,i]]
			points(jitter(array(0,dim=c(1,length(temp_data))),10)+i,
                   temp_data,col=colors[[i]], pch=pch, cex=point.cex, )
		}
	}
    if (with.p.values) {
        data_set_lengths = lapply(data,length);
        if (any(data_set_lengths > 20000)) {
            conf_data = wilcox.test(data[[1]], data[[2]], correct=FALSE);
            if (conf_data$p.value <= 2.2e-16) {
                p_val_text = 'p<2.2e-16';
            } else {
                p_val_text = paste('p=',conf_data$p.value, sep='');
            }
        } else {
            conf_data = determine_median_p_value(data[[1]], data[[2]], bootstrap.rep=bootstrap.rep);
            if (length(conf_data$p.value) == 1) {
                p_val_text = paste('p<',conf_data$p.value[1], sep='');
            } else {
                p_val_text = paste(conf_data$p.value[1],'<p<',conf_data$p.value[2], sep='');
            }
        }

        x_pos = (plot_dims[2] - plot_dims[1])*p.values.pos[1] + plot_dims[1]
        y_pos = (plot_dims[4] - plot_dims[3])*p.values.pos[2] + plot_dims[3]
        
        text(x_pos,y_pos, p_val_text, col=p.values.color);
    }
    if (return_output) {
        return(box.data);
    }
}

print_ratio_conf_string <- function(data_1,data_2) {
        p_vals = determine_median_p_value(data_1, data_2);
        
        median_ratio_low = sprintf('%.0f', 100*p_vals$ratio_conf[1]);
        median_ratio_high = sprintf('%.0f', 100*p_vals$ratio_conf[2]);
        
        median_ratio_average = (p_vals$median_vals[2] - p_vals$median_vals[1])/p_vals$median_vals[1];
        median_ratio_average = sprintf('%.0f', 100*median_ratio_average);

        print(paste(median_ratio_average, "% (", median_ratio_low, "%,", median_ratio_high, "%)", sep=''));
        return(paste(median_ratio_average, "% (", median_ratio_low, "%,", median_ratio_high, "%)", sep=''));
}

hist_with_percents <- function(data, ...) {
	hist_data = hist(data, ...);

	for (i in 1:length(hist_data$counts)) {
		y_pos = 0;

		if ((hist_data$counts[i] + 0.05*max(hist_data$counts)) > 0.5*max(hist_data$counts)) {
			y_pos = hist_data$counts[i] - 0.05*max(hist_data$counts);
		} else {
			y_pos = hist_data$counts[i] + 0.05*max(hist_data$counts);
		}
		text(hist_data$mids[i], y_pos, sprintf('%.03f',hist_data$counts[i]/sum(hist_data$counts)), srt = 45)
	}
}

get_legend_rect_points <- function(left_x,bottom_y,right_x,top_y,box_num) {
	left_x_seq = array(left_x,box_num)
	right_x_seq = array(right_x,box_num)
	bottom_y_seq = c()
	top_y_seq = c()
	
	for (i in 1:box_num) {
		bottom_y_seq = c(bottom_y_seq, (top_y - bottom_y)*((i-1)/11)+bottom_y)
		top_y_seq = c(top_y_seq,(top_y - bottom_y)*(i/11)+bottom_y)
	}
	rbind(left_x_seq,bottom_y_seq,right_x_seq,top_y_seq)
}

plot_signif_bracket <- function(upper_left, lower_right, orientation = 'regular', 
    over_text = NA, cex_text = 1.5, text_sep = 0.5, text_x_adj=0) {
    if (orientation == 'regular') {
        lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
              c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
        if (! is.na(over_text)) {
            text(mean(c(upper_left[1],lower_right[1])),upper_left[2]*1.025,over_text,cex=cex_text)
        }
    } 
    if (orientation == 'upside_down') {
        lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
              c(upper_left[2],lower_right[2],lower_right[2],upper_left[2]))  
        if (! is.na(over_text)) {
            text(mean(c(upper_left[1],lower_right[1]))+text_x_adj,lower_right[2]-text_sep,over_text,cex=cex_text)
        }
    }
}

plot_stage_length_data <- function(stage_length_data, type='stacked', top_gap = 0.1, names = NA, 
    sideways_high_xlim = 12) {

    bar_lengths = stage_length_data$bar_lengths;
    conf_ints = stage_length_data$conf_ints;

    if (all(is.na(names))) {
        names = c('Exp 1', 'Exp 2');
    }
    if (type == 'stacked') {
        conf_ints[1:3,1] = 0.22;
        conf_ints[4:6,1] = 0.57;

        err_bars = conf_ints;
        err_bars[1,3:4] = err_bars[1,3:4];
        err_bars[2,3:4] = err_bars[2,3:4] + sum(bar_lengths[1,1]);
        err_bars[3,3:4] = err_bars[3,3:4] + sum(bar_lengths[1:2,1]);
        err_bars[4,3:4] = err_bars[4,3:4];
        err_bars[5,3:4] = err_bars[5,3:4] + sum(bar_lengths[1,2]);
        err_bars[6,3:4] = err_bars[6,3:4] + sum(bar_lengths[1:2,2]);

        barplot(bar_lengths, names=names, 
            ylab='Time (min)', width=matrix(0.3,3,2), xlim=c(0,1),
            legend=c('Assembly','Stability','Disassembly'),ylim = c(0,max(err_bars)+top_gap));
        errbar(err_bars[,1],err_bars[,2],err_bars[,3],err_bars[,4],add=TRUE,cex=1E-20, xlab='', ylab='');
        return(err_bars);
    } else {
        barplot(t(bar_lengths), names=c('Assembly','Stability','Disassembly'), 
            beside=TRUE, xlim=c(0,sideways_high_xlim),
            ylab='Time (min)', legend=names,
            ylim = c(0,max(conf_ints, na.rm=T)+top_gap));

        sideways_err_bars = conf_ints;
        sideways_err_bars[1,1] = 1.5;
        sideways_err_bars[2,1] = 4.5;
        sideways_err_bars[3,1] = 7.5;
        sideways_err_bars[4,1] = 2.5;
        sideways_err_bars[5,1] = 5.5;
        sideways_err_bars[6,1] = 8.5;

        errbar(sideways_err_bars[,1],sideways_err_bars[,2],sideways_err_bars[,3],
            sideways_err_bars[,4],add=TRUE,cex=1E-20, xlab='', ylab='');
        return(sideways_err_bars);
    }
}

########################################
#Data Summary functions
########################################
filter_results <- function(results, min.r.sq=0.9, max.p.val = 0.05, debug = FALSE,
        pos.slope = TRUE) {

    points = list();
    for (i in 1:length(results)) {
        if (debug) {
            print(paste("Working on Experiment #:",i));
        }
        res = results[[i]];

        filter_sets = produce_rate_filters(res, min.r.sq = min.r.sq, max.p.val = max.p.val, pos.slope = pos.slope)
        
        assembly_filt = filter_sets$assembly;
        disassembly_filt = filter_sets$disassembly;
        joint_filt = filter_sets$joint;
        
        vars_to_add = names(res$assembly)
        for (this_var in vars_to_add) {
            points$assembly[[this_var]] = c(points$assembly[[this_var]], 
                res$assembly[[this_var]][assembly_filt]);

        }
        points$assembly$exp_num = c(points$assembly$exp_num, rep(i,length(which(assembly_filt))));
        points$assembly$lin_num = c(points$assembly$lin_num, which(assembly_filt));
        
        vars_to_add = names(res$dis)
        for (this_var in vars_to_add) {
            points$disassembly[[this_var]] = c(points$disassembly[[this_var]], 
                res$disassembly[[this_var]][disassembly_filt]);

        }
        points$disassembly$exp_num = c(points$disassembly$exp_num, rep(i,length(which(disassembly_filt))));
        points$disassembly$lin_num = c(points$disassembly$lin_num, which(disassembly_filt));
 
        vars_to_add = names(res$dis)
        for (this_var in vars_to_add) {
            points$disassembly[[this_var]] = c(points$disassembly[[this_var]], 
                res$disassembly[[this_var]][disassembly_filt]);

        }
        points$disassembly$exp_num = c(points$disassembly$exp_num, rep(i,length(which(disassembly_filt))));
        points$disassembly$lin_num = c(points$disassembly$lin_num, which(disassembly_filt));
        
        points$joint$exp_num = c(points$joint$exp_num, rep(i,length(which(joint_filt))));
        points$joint$lin_num = c(points$joint$lin_num, which(joint_filt));
        points$joint$assembly_length = c(points$joint$assembly_length, res$assembly$length[joint_filt]);
        points$joint$disassembly_length = c(points$joint$disassembly_length, res$dis$length[joint_filt]);
        
        # points$joint$assembly = c(points$joint$assembly, res$assembly$slope[joint_filt]);
        # points$joint$assembly_R = c(points$joint$assembly_R, res$assembly$R_sq[joint_filt]);
        # points$joint$disassembly = c(points$joint$disassembly, res$disassembly$slope[joint_filt]);
        # points$joint$disassembly_R = c(points$joint$disassembly_R, res$disassembly$R_sq[joint_filt]);

        # points$joint$stable_lifetime = c(points$joint$stable_lifetime, res$stable_lifetime[joint_filt]);
        # points$joint$total_lifetime = c(points$joint$total_lifetime, res$exp_props$longevity[joint_filt]);
        # points$joint$stable_mean = c(points$joint$stable_mean, res$stable_mean[joint_filt]);
        # points$joint$stable_variance = c(points$joint$stable_variance, res$stable_variance[joint_filt]);
        # for (j in names(points$joint)) {
        #     if (length(points$joint[[1]]) != length(points$joint[[j]])) {
        #         print(paste("Length of property '",j, "' wrong in joint filtering", sep=''))
        #         print(paste("Lengths are", length(points$joint[[1]]), " ", length(points$joint[[j]])))
        #         stop()
        #     }
        # }
        # 
        # if (! all(is.na(cell_intensities))) {
        #     points$control$assembly_slope = c(points$control$assembly_slope,  mean(res$assembly$slope[assembly_filt]));
        #     points$control$disassembly_slope = c(points$control$disassembly_slope,  mean(res$disassembly$slope[disassembly_filt]));
        #     points$control$mean_cell_int = c(points$control$mean_cell_int, mean(cell_intensities[[i]]));
        # }
    }
    
    points$assembly$exp_num = as.factor(points$assembly$exp_num);
    points$disassembly$exp_num = as.factor(points$disassembly$exp_num);
    points$joint$exp_num = as.factor(points$joint$exp_num);

    points$assembly = as.data.frame(points$assembly);
    points$disassembly = as.data.frame(points$disassembly);
    points$joint = as.data.frame(points$joint);

    points
}

determine_death_rate <- function(lineage_time_series) {
    total_deaths = sum(lineage_time_series$exp_props$death_status);
    total_time = dim(lineage_time_series$exp_data)[[2]];

    total_deaths/total_time;
}

determine_birth_rate <- function(lineage_time_series) {
    filtered_exp_data = lineage_time_series$exp_data[! lineage_time_series$exp_props$split_birth_status,1]

    total_births = sum(is.nan(filtered_exp_data));
    total_time = dim(lineage_time_series$exp_data)[[2]];

    total_births/total_time;
}

gather_barplot_properties <- function(data_sets, bootstrap.rep = 50000) {
    plot_props = list()
    if (is.numeric(data_sets)) {
        plot_props$mean = mean(data_sets);
        temp_conf_data = determine_mean_conf_int(data_sets, bootstrap.rep);

        plot_props$yminus = temp_conf_data[1];
        plot_props$yplus = temp_conf_data[2];
    } else {
        for (i in 1:length(data_sets)) {
            plot_props$mean = c(plot_props$mean, mean(data_sets[[i]]));

            temp_conf_data = determine_mean_conf_int(data_sets[[i]], bootstrap.rep);

            plot_props$yminus = c(plot_props$yminus, temp_conf_data[1]);
            plot_props$yplus = c(plot_props$yplus, temp_conf_data[2]);
        }
    }
    return(plot_props);
}

determine_mean_conf_int <- function(data, bootstrap.rep = 50000) {
	require(boot);
	boot_samp = boot(data, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep);
		
    conf_int = boot.ci(boot_samp, type="perc", conf=0.95)$percent[4:5]

    return(conf_int)
}

determine_mean_p_value <- function(data_1,data_2, bootstrap.rep = 50000) {
	require(boot);
	boot_samp_1 = boot(data_1, function(data_1,indexes) mean(data_1[indexes],na.rm=T), bootstrap.rep);
	boot_samp_2 = boot(data_2, function(data_1,indexes) mean(data_1[indexes],na.rm=T), bootstrap.rep);
		
    conf_int_one = boot.ci(boot_samp_1, type="perc", conf=0.99)
    conf_int_two = boot.ci(boot_samp_2, type="perc", conf=0.99)

    results = list();

    results$p_val = find_p_val_from_bootstrap(boot_samp_1, boot_samp_2);
    results$mean_vals = c(boot_samp_1$t0, boot_samp_2$t0);
    results$conf_ints = rbind(conf_int_one$percent[4:5], conf_int_two$percent[4:5])

    return(results);
}

determine_median_p_value <- function(data_1,data_2, bootstrap.rep = 50000) {
	require(boot);
    
    data_package = list(one = data_1, two = data_2);
    boot_ratio = boot(data_package, function(values, indexes) ratio_samp(values$one, values$two), bootstrap.rep);
    boot_ratio_conf = boot.ci(boot_ratio,type="perc", conf=0.99)

    results = list()
	
    results$ratio_conf = boot_ratio_conf$perc[4:5];
	
    boot_samp_1 = boot(data_1, function(values,indexes) median(values[indexes],na.rm=T), bootstrap.rep);
	boot_samp_2 = boot(data_2, function(values,indexes) median(values[indexes],na.rm=T), bootstrap.rep);
	results$p.value = find_p_val_from_bootstrap(boot_samp_1, boot_samp_2);
	results$median_vals = c(boot_samp_1$t0, boot_samp_2$t0);
    
	return(results);
}

ratio_samp <- function(data_1, data_2) {
    no_na_data_1 = na.omit(data_1);
    no_na_data_2 = na.omit(data_2);
    samp_1 = sample(no_na_data_1, length(no_na_data_1), replace=TRUE);
    samp_2= sample(no_na_data_2, length(no_na_data_2), replace=TRUE);

    return((median(samp_2) - median(samp_1))/median(samp_1));
}

gather_stage_lengths <- function(results_1, results_2, bootstrap.rep = 50000, debug=FALSE) {
	require(boot);
	bar_lengths = matrix(NA,3,2);
	conf_ints = matrix(NA,6,4);
	p_vals = rep(NA,3);
	counts = matrix(NA,3,2, dimnames=list(c('Assembly', 'Stability', 'Disassembly'), c('Results_1','Results_2')));
	
	#Assembly phase lengths
	boot_samp = boot(results_1$a$length, function(data,indexes) mean(data[indexes], na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp,type="perc")
	conf_ints[1,3:4] = boot_conf$perc[4:5]
	conf_ints[1,2] = boot_conf$t0
	bar_lengths[1,1] = boot_conf$t0
	counts[1,1] = length(results_1$a$length)
	if (debug) {
		print('Done 1');
	}
	
	boot_samp_2 = boot(results_2$a$length, function(data,indexes) mean(data[indexes], na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp_2,type="perc")
	conf_ints[4,3:4] = boot_conf$perc[4:5]
	conf_ints[4,2] = boot_conf$t0
	bar_lengths[1,2] = boot_conf$t0
	counts[1,2] = length(results_2$a$length)
	if (debug) {
		print('Done 2');
	}
	p_vals[1] = find_p_val_from_bootstrap(boot_samp,boot_samp_2)

	#Stability phase lengths
	boot_samp = boot(results_1$joint$stable_lifetime, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp,type="perc")
	conf_ints[2,3:4] = boot_conf$perc[4:5]
	conf_ints[2,2] = boot_conf$t0
	bar_lengths[2,1] = boot_conf$t0
	counts[2,1] = length(results_1$joint$stable_lifetime)
	if (debug) {
		print('Done 3');
	}

	boot_samp_2 = boot(results_2$joint$stable_lifetime, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp_2,type="perc")
	conf_ints[5,3:4] = boot_conf$perc[4:5]
	conf_ints[5,2] = boot_conf$t0
	bar_lengths[2,2] = boot_conf$t0
	counts[2,2] = length(results_2$joint$stable_lifetime)
	if (debug) {
		print('Done 4');
	}
	p_vals[2] = find_p_val_from_bootstrap(boot_samp,boot_samp_2)
	
	#Disassembly Stage Lengths
	boot_samp = boot(results_1$d$length, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp,type="perc")
	conf_ints[3,3:4] = boot_conf$perc[4:5]
	conf_ints[3,2] = boot_conf$t0
	bar_lengths[3,1] = boot_conf$t0
	counts[3,1] = length(results_1$d$length)
	if (debug) {
		print('Done 5');
	}

	boot_samp_2 = boot(results_2$d$length, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp_2,type="perc")
	conf_ints[6,3:4] = boot_conf$perc[4:5]
	conf_ints[6,2] = boot_conf$t0
	bar_lengths[3,2] = boot_conf$t0
	counts[3,2] = length(results_2$d$length)
	if (debug) {
		print('Done 6');
	}
	p_vals[3] = find_p_val_from_bootstrap(boot_samp,boot_samp_2)
		
	rm(boot_samp, boot_samp_2)
	gc()
	
    # counts = data.frame(counts, row.names = c('Assembly', 'Stability', 'Disassembly'));

	return_data <- list(bar_lengths = bar_lengths, conf_ints = conf_ints, p_vals = p_vals, counts = counts)
	return_data
}

find_p_val_from_bootstrap <- function(boot_one, boot_two, 
	p_vals_to_test = c(seq(0.99,0.01,by=-0.01),1E-3,1E-4,1E-5)) {
	
	stopifnot(class(boot_one) == "boot")
	stopifnot(class(boot_two) == "boot")
	
	correct_p_vals = c();
    for (this_val in p_vals_to_test) {
		conf_int_one = boot.ci(boot_one, type="perc", conf=(1-this_val))
		conf_int_two = boot.ci(boot_two, type="perc", conf=(1-this_val))
				
		if (ranges_overlap(conf_int_one$perc[4:5], conf_int_two$perc[4:5])) {
            break;
        } else {
            correct_p_vals = c(correct_p_vals,this_val);
		}
	}
    
    if (length(correct_p_vals) == length(p_vals_to_test)) {
        correct_p_vals = correct_p_vals[length(correct_p_vals)];
    } else {
        correct_p_vals = correct_p_vals[(length(correct_p_vals)-1):length(correct_p_vals)];
    }

    return(correct_p_vals)
}

ranges_overlap <- function(range_1, range_2) {
	stopifnot(length(range_1) == 2)
	stopifnot(length(range_2) == 2)
	
	if (range_1[1] < range_2[1] & range_1[2] > range_2[1]) {
		return(TRUE);
	} else if (range_1[1] > range_2[1] & range_1[1] < range_2[2]) {
		return(TRUE);
	}
	return(FALSE);
}

gather_general_dynamic_props <- function(results, debug=FALSE) {
	points = list()
	for (i in 1:length(results)) {
		res = results[[i]]
        if (debug) {
            print(paste("Working on", i))
        }

        filt = ! res$exp_props$split_birth_status & res$exp_props$death_status

		points$longevity = c(points$longevity, res$exp_props$longevity[filt])
		points$mean_area = c(points$mean_area, res$exp_props$mean_area[filt])
		points$mean_axial_ratio = c(points$mean_axial_ratio, res$exp_props$mean_axial_ratio[filt])
		points$mean_edge_dist = c(points$mean_edge_dist, res$exp_props$mean_edge_dist[filt])
		points$largest_area = c(points$largest_area, res$exp_props$largest_area[filt])
		points$ad_sig = c(points$ad_sig, res$exp_props$ad_sig[filt])
		points$average_speed = c(points$average_speed, res$exp_props$average_speeds[filt])
		points$exp_num = c(points$exp_num, rep(i,length(which(filt))))
	}
    
    points$exp_num = as.factor(points$exp_num);

	points = as.data.frame(points)
	points
}

gather_static_props <- function(ind_results, debug=FALSE) {
	ind_data = list();
	
	for (i in 1:length(ind_results)) {
		res = ind_results[[i]]
        if (debug) {
            print(paste("Working on", i))
        }
		filt_by_area = res$Area >= min(res$Area)# & res$I_num == 1
		ind_data$Area = c(ind_data$Area, res$Area[filt_by_area]);
		ind_data$ad_sig = c(ind_data$ad_sig, res$Average_adhesion_signal[filt_by_area]);
		ind_data$ad_var = c(ind_data$ad_var, res$Variance_adhesion_signal[filt_by_area]);
		ind_data$axial_r = c(ind_data$axial_r, res$MajorAxisLength[filt_by_area]/res$MinorAxisLength[filt_by_area]);
	
		ind_data$cent_dist = c(ind_data$cent_dist, res$Centroid_dist_from_edge[filt_by_area]);
        ind_data$exp_num = c(ind_data$exp_num, rep(i, length(filt_by_area)));
	}
    
    ind_data$exp_num = as.factor(ind_data$exp_num);
	as.data.frame(ind_data)
}

################################################################################
# File Reading/Writing Functions
################################################################################
load_results <- function(dirs,file, debug=TRUE) {
	results = list()
    if (length(dirs) == 0) return()

	for (i in 1:length(dirs)) {
		this_file = file.path(dirs[i],file)
		if (file.exists(this_file)) {
            if (debug) {
                print(paste('Loading File:', this_file))
            }
			var_name = load(file.path(dirs[i],file))
			results[[i]] = get(var_name)
		}
	}
	results
}

load_data_files <- function(dirs, files, headers = FALSE, inc_exp_names = TRUE, debug = FALSE) {
    #This function searches through the directories provided "dirs" and looks
    #for the files "files". If all the files are present in all the dirs, they
    #are loaded. The headers are included if the matching headers position
    #is true, otherwise, the headers are excluded
    results = list()
    if (length(dirs) == 0) return()
	
    exp_names = list()
	all_files_present_dirs = c()	
	for (i in 1:length(dirs)) {
		seen_files = 0
		all_files_present = 1;
		for (j in 1:length(files)) {
			this_file = file.path(dirs[i],files[j])
		
			if (! file.exists(this_file)) {
				all_files_present = 0;
			}
		}
		if (all_files_present) {
			all_files_present_dirs = c(all_files_present_dirs, dirs[i]);
			
			regex_range = regexpr("time_series_[[:digit:]]",dirs[[i]])
	        if (regex_range[1] == -1) {
        		exp_names[[length(all_files_present_dirs)]] = dirs[[i]];
	        } else {
    	    	exp_names[[length(all_files_present_dirs)]] = substr(dirs[[i]], regex_range[1], regex_range[1] + attr(regex_range,'match.length'));
        	}
		}
	}
	
	if (debug) {
        print("Directories with full sets of files:");
		print(all_files_present_dirs)
        print('');
	}

    stopifnot(length(all_files_present_dirs) != 0)
        
	
	for (i in 1:length(files)) {
		results[[i]] = list()
	}
	
	header_array = array(headers, dim=c(1,length(all_files_present_dirs)))

	for (i in 1:length(all_files_present_dirs)) {
        if (debug) {
            print(paste("Working on directory:",all_files_present_dirs[i]))
        }
		for (j in 1:length(files)) {
            if (debug) {
                print(paste("Working on file:",files[j]))
                print(j)
            }
		    results[[j]][[i]] = read.table(file.path(all_files_present_dirs[i],files[j]), header = header_array[i], sep=",")
            # if (headers) {
            # } else {
			    # results[[j]][[i]] = read.table(file.path(all_files_present_dirs[i],files[j]), header = header_array[i], sep=",")
            # }
		}
	}
	if (inc_exp_names) {
		results[[length(files) + 1]] = exp_names
	}
		
	if (length(results) == 1) {
		results = results[[1]]
	}
	results
}

########################################
#Spacial Functions
########################################
correlate_signal_vs_dist <- function(intensities, cent_x, cent_y, min_overlap=40) {
	correlations = c()
	distances = c()
	for (i in 1:dim(intensities)[[1]]) {		
		this_row = intensities[i,];
		if (sum(! is.na(this_row)) < min_overlap) {
			next
		}
		overlap_rows = which(apply(intensities,1, function(data_2) enough_overlap(this_row, data_2, min_overlap=min_overlap)));		
		overlap_rows = setdiff(overlap_rows,1:i)
		
		for (j in overlap_rows) {
			overlap_entries = which(! is.na(this_row) & ! is.na(intensities[j,]))
			correlations = c(correlations, cor(as.numeric(this_row[overlap_entries]), as.numeric(intensities[j,overlap_entries])))
			
			temp_dists = c()
			cent_x_start = cent_x[i,]
			cent_y_start = cent_y[i,]
			cent_x_end = cent_x[j,]
			cent_y_end = cent_y[j,]
			for (k in overlap_entries) {
				temp_dists = c(temp_dists, sqrt((cent_x_start[k] - cent_x_end[k])^2 + (cent_y_start[k] - cent_y_end[k])^2))
			}
			distances = c(distances, mean(temp_dists))
		}
	}
	return(list(correlations = correlations, distances = distances))
}

enough_overlap <- function(data_1, data_2, min_overlap=10) {
	if (length(data_1) != length(data_2)) {
		return(FALSE)
	}
	
	if (length(which(! is.na(data_1) & ! is.na(data_2))) >= min_overlap) {
		return(TRUE);
	} else {
		return(FALSE);
	}	
}
stopifnot(! enough_overlap(rep(NaN, 10), rep(1,10), min_overlap=10))
stopifnot(enough_overlap(rep(2, 10), rep(1,10), min_overlap=10))
stopifnot(! enough_overlap(c(rep(1, 9), NaN), rep(1,10), min_overlap=10))

bin_corr_data <- function(corr_results, bin_size = NA, bootstrap.rep = 5000, bin_max = NA, pixel_size=0.215051) {
	require(boot)
	corr_results$distances = corr_results$distances*pixel_size

	if (is.na(bin_size)) {
		bin_size = pixel_size * 5
	}

	if (! is.na(bin_max)) {
		bins = c(0,seq(bin_size, bin_max*pixel_size, by=bin_size))
	} else {
		bins = c(0,seq(bin_size,max(corr_results$distances)+bin_size, by=bin_size))
	}
	means = rep(NA, length(bins) - 1);
	upper = rep(NA, length(bins) - 1);
	lower = rep(NA, length(bins) - 1);
	counts = rep(NA, length(bins) - 1);
	bin_mids = rep(NA, length(bins) - 1);
	
	for (i in 2:length(bins)) {
		this_bin_data = corr_results$correlations[corr_results$distances > bins[i-1] & corr_results$distances <= bins[i]]
		means[i-1] = mean(this_bin_data)
		bin_mids[i-1] = mean(c(bins[i - 1], bins[i]))
		counts[i-1] = length(this_bin_data)
		if (counts[i-1] < 10) {
			next
		}
		boot_samp = boot(this_bin_data, function(data,indexes) mean(data[indexes]), bootstrap.rep)
		if (all(boot_samp$t[1] == boot_samp$t)) {
			browser()
		} else {		
			boot_conf = boot.ci(boot_samp,type="perc")
			lower[i-1] = boot_conf$perc[4]
			upper[i-1] = boot_conf$perc[5]
		}
	}
	
	return(as.data.frame(list(means = means, mids = bin_mids, upper = upper, lower = lower, counts = counts)))
}

################################################################################
# Testing
################################################################################
stopifnot(! ranges_overlap(c(0,1),c(-0.5,-0.01)))
stopifnot(ranges_overlap(c(0,1),c(-0.5,0.01)))
stopifnot(ranges_overlap(c(0,1),c(0.5,0.9)))
stopifnot(! ranges_overlap(c(0,1),c(1.01,1.2)))

stopifnot(! ranges_overlap(c(-0.5,-0.01),c(0,1)))
stopifnot(ranges_overlap(c(-0.5,0.01),c(0,1)))
stopifnot(ranges_overlap(c(0.5,0.9),c(0,1)))
stopifnot(! ranges_overlap(c(1.01,1.2),c(0,1)))

################################################################################
# Main Program
################################################################################

args = commandArgs(TRUE);
if (length(args) != 0) {
    debug = FALSE;
    min_length = 10;
    
	#split out the arguments from the passed in parameters and assign variables 
	#in the current scope
    for (this_arg in commandArgs()) {
        split_arg = strsplit(this_arg,"=",fixed=TRUE)
            if (length(split_arg[[1]]) == 1) {
                assign(split_arg[[1]][1], TRUE);
            } else {
                assign(split_arg[[1]][1], split_arg[[1]][2]);
            }
    }
	
    class(min_length) <- "numeric";
	print(data_dir);
    if (exists('data_dir') & exists('model_type')) {
        if (model_type == 'average') {
            average_model = gather_bilinear_models_from_dirs(data_dir,
                    data_file='Average_adhesion_signal.csv', min_length = min_length,
                    results.file=file.path('..','models','intensity.Rdata'), debug=debug)
            print(class(average_model))
            write_assembly_disassembly_periods(average_model[[1]],file.path(data_dir,'..'))	
        }
        if (model_type == 'cell_background') {
            temp = gather_bilinear_models_from_dirs(data_dir, 
                    data_file='CB_corrected_signal.csv', min_length = min_length,
                    results.file=file.path('..','models','CB_corrected.Rdata'), debug=debug)
        }
        if (model_type == 'local_background') {
            temp = gather_bilinear_models_from_dirs(data_dir, 
                    data_file='Background_corrected_signal.csv', min_length = min_length,
                    results.file=file.path('..','models','local_corrected.Rdata'), debug=debug)
        }
        if (model_type == 'area') {
            temp = gather_bilinear_models_from_dirs(data_dir, 
                    data_file='Area.csv', min_length = min_length, log.trans = FALSE,
                    results.file=file.path('..','models','area.Rdata'), debug=debug)
        }
        if (model_type == 'box_intensity') {
            temp = gather_bilinear_models_from_dirs(data_dir, 
                    data_file='Box_intensity.csv', min_length = min_length,
                    results.file=file.path('..','models','box.Rdata'), debug=TRUE)
        }
        if (model_type == 'background_correlation_model') {
            intensity_data <- 
                read.table(file.path(data_dir,'Background_corrected_signal.csv'), 
                        sep=',', header=FALSE);

            centroid_x <- read.table(file.path(data_dir,'Centroid_x.csv'), sep=',', header=FALSE);
            centroid_y <- read.table(file.path(data_dir,'Centroid_y.csv'), sep=',', header=FALSE);

            results = correlate_signal_vs_dist(intensity_data, centroid_x, centroid_y)
                output_file = file.path(data_dir,'..','models','background_corr.Rdata')
                if (! file.exists(dirname(output_file))) {
                    dir.create(dirname(output_file),recursive=TRUE)
                }
            save(results,file = output_file);
        }
    }
}
