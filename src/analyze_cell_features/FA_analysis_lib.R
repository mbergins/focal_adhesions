################################################################################
# FA_analysis_lib.R: various functions used to work with the focal adhesion data
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
		
		plot(0:(length(ad_seq)-1), ad_seq, xlab='Time (minutes)', ylab='Adhesion Intensity',type="o");
		
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

plot_ad_intensity <- function (results,index,phase_lengths,R_sq,y_limits,slope,time.spacing=1, ...) {
    
    #pull out the intensity signal of interest and build the time vector for
    #plotting
	ad_seq = as.vector(results$exp_data[index,])
    time_points = which(!is.nan(ad_seq))*time.spacing - time.spacing
    ad_seq = na.omit(ad_seq)
    
    plot(time_points, ad_seq, xlab='Time (minutes)', ylab='Average Intensity',type="o", 
        main=index,ylim=y_limits);
    
    #add bars at the first point of the movie and the last, to indicate why a
    #slope value wasn't calculated
    plot_limits = par("usr")
    segments(0,par("usr")[3],0,par("usr")[4])
    segments((length(results$exp_data[index,])-1)*time.spacing,par("usr")[3],
        (length(results$exp_data[index,])-1)*time.spacing,par("usr")[4])
    
    if (is.finite(phase_lengths[1])) {
        assemb_time = time_points[1:phase_lengths[1]]
        assemb_sig = ad_seq[1:phase_lengths[1]]
        log_model = lm(log(assemb_sig) ~ assemb_time);
        predicted_sig = predict(log_model);
        
        lines(assemb_time, exp(predicted_sig), col='darkgreen',lwd = 3)
        if (is.finite(R_sq[1])) {
            text(plot_limits[1],plot_limits[4]*0.95,pos=4, 
                paste('R-sq:', sprintf('%.2f',R_sq[1])))
        }
        if (is.finite(slope[1])) {
            text(plot_limits[1],plot_limits[4]*0.925,pos=4, 
                paste('Rate:', sprintf('%.4f',slope[1])))
        }
    }
    
    if (is.finite(phase_lengths[2])) {
        dis_time = rev(time_points)[1:phase_lengths[2]]
        dis_sig = rev(ad_seq)[1:phase_lengths[2]] 
        log_model = lm(log(dis_sig) ~ dis_time);
        predicted_sig = predict(log_model);
        lines(dis_time, exp(predicted_sig), col='red',lwd = 3)
        
        if (is.finite(R_sq[2])) {
            text(plot_limits[2],plot_limits[4]*0.95,pos=2, 
                paste('R-sq:', sprintf('%.2f',R_sq[2])))
        }
        if (is.finite(slope[2])) {
            text(plot_limits[2],plot_limits[4]*0.925,pos=2, 
                paste('Rate:', sprintf('%.4f',slope[2])))
        }
    }
}

plot_all_ad_intensities <- function (raw_data,filtered_set) {
    exp_count = length(raw_data)

    for (i in 1:exp_count) {
        this_raw_data = raw_data[[i]];
        this_assembly_set = subset(filtered_set$assembly, exp_num == i);
        this_disass_set = subset(filtered_set$dis, exp_num == i);

        sorted_uniq_lin_nums = sort(unique(c(this_disass_set$lin_num,this_assembly_set$lin_num)))
        
        pdf(file.path(this_raw_data$exp_dir,'ad_intensity_plots.pdf'))
        for (ad_num in sorted_uniq_lin_nums) {
            as_row = which(this_assembly_set$lin_num == ad_num)
            dis_row = which(this_disass_set$lin_num == ad_num)
            
            phase_lengths = c(NA, NA)
            R_sq = c(NA, NA)
            if (length(as_row) > 0) {
                phase_lengths[1] = this_assembly_set$length[as_row]
                R_sq[1] = this_assembly_set$adj.r.squared[as_row]
            }
            if (length(dis_row) > 0) {
                phase_lengths[2] = this_disass_set$length[dis_row]
                R_sq[2] = this_disass_set$adj.r.squared[dis_row]
            }
            
            plot_ad_intensity(this_raw_data,ad_num,phase_lengths,R_sq,time.spacing = 2.5)
        }
        graphics.off()
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
    colors=c('darkgreen','red','brown','blue','yellow','pink','cyan','gray','orange','purple'),
    notch=F,  range=1.5, inc.n.counts = TRUE, inc.points = TRUE, pch=20,
    na.omit = TRUE, point.cex=0.2, return_output = FALSE, with.p.value = TRUE,
    p.value.pos = c(0.5,0.9), p.value.color = 'black', p.value.type = 'mean', 
    bootstrap.rep = 10000, ...) {
	
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
	
	box.data = boxplot(data,notch = notch,names = names,range = range,pch=19,cex=0.25,...)

    plot_dims = par("usr");
	if (inc.points) {
        colors = rep(colors,length(data));
		for (i in 1:length(data)) {
			this_data = data[[i]]
			temp_data = this_data[this_data >= box.data$stat[1,i] & this_data <= box.data$stat[5,i]]
			points(jitter(array(0,dim=c(1,length(temp_data))),10)+i,
                   temp_data,col=colors[i], pch=pch, cex=point.cex)
		}
	}
    if (with.p.value) {
        data_set_lengths = lapply(data,length);
        if (any(data_set_lengths > 20000)) {
            conf_data = wilcox.test(data[[1]], data[[2]], correct=FALSE);
            if (conf_data$p.value <= 2.2e-16) {
                p_val_text = 'p<2.2e-16';
            } else {
                p_val_text = paste('p=',conf_data$p.value, sep='');
            }
        } else {
            if (p.value.type == 'median') {
                conf_data = determine_median_p_value(data[[1]], data[[2]],bootstrap.rep);
            } else if (p.value.type == 'mean') {
                conf_data = determine_mean_p_value(data[[1]], data[[2]],bootstrap.rep);
            } else {
                print('Unrecognized p.value.type');
                stop()
            }
            p_val_text = paste('p<',conf_data$p.value[1], sep='');
            if (conf_data$p.value[1] < 0.05) {
                if (p.value.type == 'median') {
                    p_val_text = paste(p_val_text,'\n',
                        round(100*(median(data[[2]])/median(data[[1]]) - 1)), '%',sep='')
                } else {
                    p_val_text = paste(p_val_text,'\n',
                        round(100*(mean(data[[2]])/mean(data[[1]]) - 1)), '%',sep='')
                }
            }
        }

        x_pos = (plot_dims[2] - plot_dims[1])*p.value.pos[1] + plot_dims[1]
        y_pos = (plot_dims[4] - plot_dims[3])*p.value.pos[2] + plot_dims[3]
        
        text(x_pos,y_pos, p_val_text, col=p.value.color);
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
    over_text = NA, cex_text = 1,text_sep = 0, text_x_adj=0) {
    if (orientation == 'regular') {
        lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
              c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
        if (! is.na(over_text)) {
            text(mean(c(upper_left[1],lower_right[1])),upper_left[2]-text_sep,over_text,pos=3,cex=cex_text)
        }
    } 
    if (orientation == 'upside_down') {
        lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
              c(upper_left[2],lower_right[2],lower_right[2],upper_left[2]))  
        if (! is.na(over_text)) {
            text(mean(c(upper_left[1],lower_right[1]))+text_x_adj,lower_right[2]-text_sep,over_text,cex=cex_text)
            # text(mean(c(upper_left[1],lower_right[1])),lower_right[2],over_text,pos=3,cex=cex_text)
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

add_labels_with_sub <- function(data_sets,names,at,subtitle=NA,with.axis=T,...) {
    pl_size = par("usr");
    char_size = par("cxy")[2]
    
    stopifnot(length(at) == length(names))
    
    if (any(is.na(data_sets))) {
        labs = names;
    } else {
        labs = c()
        for (i in 1:length(data_sets)) {
            labs = c(labs, paste(names[i],' (n=',length(data_sets[[i]]),')',sep=''))
        }
    }
    
    if (with.axis) {
        axis(1,labels=labs,at=at,...)
    }

    if (! is.na(subtitle)) {
        lines(c(min(at),max(at)),rep(pl_size[3]-char_size*1.675,2),lwd=3)
        mtext(subtitle,1,at=mean(at),line=1.7)
    }
}

plot_rates_vs_time <- function(before, after) {
    plot(slope ~ birth_i_num, data=before$assembly,main='Assembly')
    abline(lm(slope ~ birth_i_num, data=before$assembly),col='red',lwd=2)

    plot(slope ~ birth_i_num, data=after$assembly,main='Assembly')
    abline(lm(slope ~ birth_i_num, data=after$assembly),col='red',lwd=2)

    plot(slope ~ death_i_num, data=before$dis,main='Disassembly')
    abline(lm(slope ~ death_i_num, data=before$dis),col='red',lwd=2)

    plot(slope ~ death_i_num, data=after$dis,main='Disassembly')
    abline(lm(slope ~ death_i_num, data=after$dis),col='red',lwd=2)
}

plot_dual_hist <- function(data_set_1,data_set_2,
    colors=c(rgb(0,0,1,0.6),rgb(1,1,0,0.6)),break.int=NA,
    xlab=NA,ylab=NA,axes=F) {

    hist_data_1 = hist(data_set_1,plot=F);
    hist_data_2 = hist(data_set_2,plot=F);
    
    if (max(hist_data_1$breaks) > max(hist_data_2$breaks)) {
        these_breaks = hist_data_1$breaks;
    } else {
        these_breaks = hist_data_2$breaks;
    }

    y_max = max(hist_data_1$counts/sum(hist_data_1$counts),hist_data_2$counts/sum(hist_data_2$counts))

    if (!is.na(break.int)) {
        these_breaks = seq(0,max(c(data_set_1,data_set_2)),by=break.int)
        if (max(these_breaks) != max(c(data_set_1,data_set_2))) {
            these_breaks = c(these_breaks, max(these_breaks)+break.int)
        }
    }
    
    hist(data_set_1,col=colors[1],breaks=these_breaks,freq=F,axes=axes,main='',
        xlab=xlab,ylab=ylab,ylim=c(0,y_max));
    hist(data_set_2,col=colors[2],add=T,breaks=these_breaks,freq=F,axes=F);
}

plot_with_shading <- function(x,y,upper_conf,lower_conf,col,...) {
    col_light = c();
    for (i in 1:dim(y)[1]) {
        col_light = c(col_light,rgb(t(col2rgb(col[i]))/255,alpha=0.5));
    }
    
    for (i in 1:dim(y)[1]) {
        if (i == 1) {
            plot(x,y[i,],col=col[i],...);
        } else {
            lines(x,y[i,],col=col[i],...);
        }
        polygon(c(x,rev(x)), c(lower_conf[i,],rev(upper_conf[i,])),col=col_light[i],border=NA)
    }
    
    #re-draw each of the primary lines so that they stand out against the
    #transparency
    for (i in 1:dim(y)[1]) {
        lines(x,y[i,],col=col[i],...);
    }
}

###########################################################
#Data Summary functions
###########################################################

#######################################
# Bootstraping
#######################################
determine_mean_conf_int <- function(data, bootstrap.rep = 10000) {
	require(boot);
	
    boot_samp = boot(data, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep);
    conf_int = boot.ci(boot_samp, type="bca", conf=0.95)$bca[4:5]
    
    boot_samp <- NULL

    return(conf_int)
}

determine_median_conf_int <- function(data, bootstrap.rep = 10000) {
	require(boot);
	
    boot_samp = boot(data, function(data,indexes) median(data[indexes],na.rm=T), bootstrap.rep);
    browser()
    conf_int = boot.ci(boot_samp, type="bca", conf=0.95)$bca[4:5]

    return(conf_int)
}

determine_mean_p_value <- function(data_1,data_2, bootstrap.rep = 10000) {
	require(boot);
    
	boot_samp_1 = boot(data_1, function(values,indexes) mean(values[indexes],na.rm=T), bootstrap.rep);
	boot_samp_2 = boot(data_2, function(values,indexes) mean(values[indexes],na.rm=T), bootstrap.rep);
	
    results = list();

	results$p.value = find_p_val_from_bootstrap(boot_samp_1, boot_samp_2);
    results$mean_vals = c(boot_samp_1$t0, boot_samp_2$t0);

    return(results);
}

determine_median_ratio_conf <- function(data_1,data_2, bootstrap.rep=10000) {
	require(boot);
    
    data_package = list(one = data_1, two = data_2);
    boot_ratio = boot(data_package, function(values, indexes) ratio_samp(values$one, values$two), bootstrap.rep);
    boot_ratio_conf = boot.ci(boot_ratio,type="bca", conf=0.99)

    return(boot_ratio_conf)
}

determine_median_p_value <- function(data_1,data_2, bootstrap.rep = 10000) {
	require(boot);
    
    results = list()
	
    boot_samp_1 = boot(data_1, function(values,indexes) median(values[indexes],na.rm=T), bootstrap.rep);
	boot_samp_2 = boot(data_2, function(values,indexes) median(values[indexes],na.rm=T), bootstrap.rep);

	results$p.value = find_p_val_from_bootstrap(boot_samp_1, boot_samp_2);
	results$median_vals = c(boot_samp_1$t0, boot_samp_2$t0);
    
	return(results);
}

find_p_val_from_bootstrap <- function(boot_one, boot_two) {
    p_vals_to_test = c(seq(0.99,0.01,by=-0.01),1E-3,1E-4,1E-5)

	stopifnot(class(boot_one) == "boot")
	stopifnot(class(boot_two) == "boot")
	
    overlap = rep(NA,length(p_vals_to_test));
    overlap[length(overlap)] = bootstrap_overlap(boot_one, boot_two, 1-p_vals_to_test[length(p_vals_to_test)])
    #if the lowest p-value doesn't produce overlapping intervals, return the
    #lowest p-value from the function, otherwise, begin binary search for
    #p-value cutoff
    if (! overlap[length(overlap)]) {
        return(p_vals_to_test[length(p_vals_to_test)])
    }
    overlap[1] = bootstrap_overlap(boot_one, boot_two, 1-p_vals_to_test[1])
    
    upper = length(p_vals_to_test)
    lower = 1
    middle = floor(mean(c(upper,lower)))

    while (middle != lower && middle != upper) {
        overlap[middle] = bootstrap_overlap(boot_one, boot_two, 1-p_vals_to_test[middle])
        if (overlap[middle]) {
            upper = middle;
            middle = floor(mean(c(upper,lower)))
        } else {
            lower = middle;
            middle = floor(mean(c(upper,lower)))
        }
    }

    overlap_indexes = which(overlap)
    no_overlap_indexes = which(!overlap)
    
    p_val_range = c(p_vals_to_test[no_overlap_indexes[length(no_overlap_indexes)]],
        p_vals_to_test[overlap_indexes[1]])
    
    # these lines deal with a strange bug where R decides that some of the
    # p-values are not quantized property, i.e. 0.05 as 0.04999999..., adding
    # this sprintf to deal with those problems
    if (p_val_range[1] > 1E-3) { p_val_range[1] = sprintf('%0.2f',as.numeric(p_val_range[1])) }
    if (p_val_range[2] > 1E-3) { p_val_range[2] = sprintf('%0.2f',as.numeric(p_val_range[2])) }
    
    return(p_val_range)
}

bootstrap_overlap <- function(boot_one,boot_two,p_val, type="bca") {
    conf_int_one = boot.ci(boot_one, type=type, conf=p_val)
    conf_int_two = boot.ci(boot_two, type=type, conf=p_val)
    
    #occassionally, very large confidence intervals yield NA values, we will
    #assume that those confidence intervals overlap
    if (any(is.na(c(conf_int_one$bca[4:5], conf_int_two$bca[4:5])))) {
        return(TRUE)
    }

    return(ranges_overlap(conf_int_one$bca[4:5], conf_int_two$bca[4:5]))
}

#######################################
# Filtering
#######################################
filter_results <- function(results, model_count = NA, min.r.sq=0.9, max.p.val = 0.05, debug = FALSE,
        pos.slope = TRUE,old.names=F,min.cent.dist=NA) {

    ad_data = list();
    
    for (i in 1:length(results)) {
        if (debug) {
            print(paste("Working on Experiment #:",i,'/',length(results)));
        }
        res = results[[i]];
        
        filter_sets = produce_rate_filters(res, model_count, min.r.sq = min.r.sq, max.p.val = max.p.val, 
            pos.slope = pos.slope,old.names=old.names,min.cent.dist=min.cent.dist)
        
        assembly_filt = filter_sets$assembly;
        disassembly_filt = filter_sets$disassembly;
        joint_filt = filter_sets$joint;
        
        prop_set = c('largest_area','ad_sig','mean_axial_ratio','birth_i_num',
            'death_i_num','mean_area','longevity','Mean_FA_cent_dist')
        if (any(names(res$exp_props) == 'drug_addition_time')) {
            prop_set = c(prop_set, 'drug_addition_time');
        }
        these_props = subset(res$exp_props, select = prop_set)
        #######################################################################
        # Building the Assembly Data Set
        #######################################################################

        this_assem_data = list()
        this_assem_data = res$assembly[assembly_filt,]
        # browser()
        this_assem_data = cbind(this_assem_data,these_props[assembly_filt,])
        this_assem_data$exp_num = c(this_assem_data$exp_num, rep(i,length(which(assembly_filt))));
        this_assem_data$lin_num = c(this_assem_data$lin_num, which(assembly_filt));

        ad_data$assembly = rbind(ad_data$assembly,this_assem_data)
             
        #######################################################################
        # Building the Disassembly Data Set
        #######################################################################
        
        this_dis_data = list()
        this_dis_data = res$disassembly[disassembly_filt,]
        this_dis_data = cbind(this_dis_data,these_props[disassembly_filt,])
        this_dis_data$exp_num = c(this_dis_data$exp_num, rep(i,length(which(disassembly_filt))));
        this_dis_data$lin_num = c(this_dis_data$lin_num, which(disassembly_filt));
        
        ad_data$disassembly = rbind(ad_data$disassembly,this_dis_data)
        
        #######################################################################
        # Building the Joint Data Set
        #######################################################################
        ad_data$joint$exp_num = c(ad_data$joint$exp_num, rep(i,length(which(joint_filt))));
        ad_data$joint$lin_num = c(ad_data$joint$lin_num, which(joint_filt));
        ad_data$joint$assembly_length = c(ad_data$joint$assembly_length, res$assem$length[joint_filt]);
        ad_data$joint$assembly_slope = c(ad_data$joint$assembly_slope, res$assem$slope[joint_filt]);
        ad_data$joint$disassembly_length = c(ad_data$joint$disassembly_length, res$dis$length[joint_filt]);
        ad_data$joint$disassembly_slope = c(ad_data$joint$disassembly_slope, res$dis$slope[joint_filt]);
        ad_data$joint$longevity = c(ad_data$joint$longevity, res$exp_props$longevity[joint_filt]);

        #the assembly/disassembly lengths are saved as number of images, not number of minutes, so to find the stability length we need the number of 
        i_num_joint_longev = res$exp_props$death_i_num[joint_filt] - res$exp_props$birth_i_num[joint_filt] + 1;
        ad_data$joint$stability_length = c(ad_data$joint$stability_length, 
            i_num_joint_longev - res$assem$length[joint_filt] - res$dis$length[joint_filt]);
    }
    
    ad_data$assembly$exp_num = as.factor(ad_data$assembly$exp_num);
    ad_data$disassembly$exp_num = as.factor(ad_data$disassembly$exp_num);
    ad_data$joint$exp_num = as.factor(ad_data$joint$exp_num);

    ad_data$assembly = as.data.frame(ad_data$assembly);
    ad_data$disassembly = as.data.frame(ad_data$disassembly);
    ad_data$joint = as.data.frame(ad_data$joint);

    ad_data
}

filter_results_by_drug_addition <- function(data_set) {
    # assembly_filter = data_set$as$birth_i_num <= data_set$as$drug_addition_time &
    #     data_set$as$death_i_num >= data_set$as$drug_addition_time
    # assembly_filter[is.na(assembly_filter)] = FALSE

    # disassembly_filter = data_set$di$birth_i_num <= data_set$di$drug_addition_time &
    #     data_set$di$death_i_num >= data_set$di$drug_addition_time
    # disassembly_filter[is.na(disassembly_filter)] = FALSE
    assembly_filter = data_set$assembly$birth_i_num < data_set$assembly$drug_addition_time &
        (data_set$assembly$birth_i_num + data_set$assembly$length) > data_set$assembly$drug_addition_time
    
    disassembly_filter = data_set$disassembly$death_i_num > data_set$disassembly$drug_addition_time &
        (data_set$d$death_i_num - data_set$d$length) < data_set$disassembly$drug_addition_time
    
    assembly_close = data_set$assembly$birth_i_num > data_set$assembly$drug_addition_time &
        data_set$assembly$birth_i_num <= data_set$assembly$drug_addition_time + 5
  
    assem_length = length(assembly_filter);
    assem_filt_length = sum(assembly_filter);
    disassem_length = length(disassembly_filter);
    disassem_filt_length = sum(disassembly_filter);

    percent_as_removed = sprintf('%.0f',assem_filt_length/assem_length*100)
    percent_di_removed = sprintf('%.0f',disassem_filt_length/disassem_length*100)
    print(paste('Removing', assem_filt_length, 
        'assembly (',percent_as_removed,'%) and', length(which(disassembly_filter)), 
        'disassembly (',percent_di_removed,'%) models'))

    data_set$assembly = data_set$assembly[!assembly_filter,]
    data_set$disassembly = data_set$disassembly[!disassembly_filter,]
    return(data_set)
}

split_results_by_drug_addition <- function(data_set) {
    assembly_filter = data_set$assembly$birth_i_num < data_set$assembly$drug_addition_time;
    
    disassembly_filter = data_set$disassembly$death_i_num < data_set$disassembly$drug_addition_time;

    final_set = list()
    final_set$assembly$before = data_set$assembly[assembly_filter,];
    final_set$assembly$after = data_set$assembly[!assembly_filter,];
    final_set$disassembly$before = data_set$disassembly[disassembly_filter,];
    final_set$disassembly$after = data_set$disassembly[!disassembly_filter,];

    return(final_set)
}

split_results_by_exp <- function(data_set) {
    per_exp = list()
    
    exp_num = union(unique(data_set$assembly$before$exp_num),unique(data_set$assembly$after$exp_num))
    temp = c() 
    for (this_exp_num in exp_num) {
        before = subset(data_set$assembly$before, exp_num==this_exp_num);
        after = subset(data_set$assembly$after, exp_num==this_exp_num);
        
        temp = c(temp, median(after$slope) - median(before$slope))
    }
    per_exp$assembly = temp;
    
    exp_num = union(unique(data_set$disassembly$before$exp_num),unique(data_set$disassembly$after$exp_num))
    temp = c() 
    for (this_exp_num in exp_num) {
        before = subset(data_set$disassembly$before, exp_num==this_exp_num);
        after = subset(data_set$disassembly$after, exp_num==this_exp_num);

        temp = c(temp, median(after$slope) - median(before$slope))
    }
    per_exp$disassembly = temp;

    return(per_exp);
}

#######################################
# Dynamics Summary
#######################################
gather_dynamics_summary <- function(data,debug=FALSE,table.text.format='%.1f\u00B1%.2f') {
    data_summary = list()

    data_summary$counts = lapply(data,length);
    data_summary$mean = lapply(data,function(x) mean(x,na.rm=T));

    data_summary$t_conf = lapply(data, function(x) t.test(x)$conf.int[1:2])
    data_summary$t_margin = lapply(data, function(x) t.test(x)$conf.int[2]-mean(x,na.rm=T))
    
    data_summary$table_text = lapply(data, 
        function(x) sprintf(table.text.format,mean(x,na.rm=T),t.test(x)$conf.int[2]-mean(x,na.rm=T)))
    
    return(data_summary);
}

gather_dynamics_summary_median <- function(data,debug=FALSE,table.text.format='%.1f') {
    data_summary = list()

    data_summary$counts = lapply(data,length);
    data_summary$mean = lapply(data,function(x) median(x,na.rm=T));

    data_summary$t_conf = lapply(data, function(x) t.test(x)$conf.int[1:2])
    data_summary$t_margin = lapply(data, function(x) t.test(x)$conf.int[2]-mean(x,na.rm=T))
    
    data_summary$table_text = lapply(data, 
        function(x) sprintf(table.text.format,median(x,na.rm=T)))
    
    return(data_summary);
}
gather_global_exp_summary <- function(data_set) {
    stopifnot(is.list(data_set))
    
    average_adhesion_count = c()
    sd_adhesion_count = c()
    exp_dirs = c()
    for (i in 1:length(data_set)) {
        these_counts = apply(data_set[[i]]$exp_data,2,function(x) sum(!is.nan(x)))
        average_adhesion_count = c(average_adhesion_count, mean(these_counts))
        sd_adhesion_count = c(sd_adhesion_count, sd(these_counts))

        exp_dirs = c(exp_dirs, data_set[[i]]$exp_dir)
    }
    
    adhesion_count_data = data.frame(mean_ad_count = average_adhesion_count, 
        sd_adhesion_count = sd_adhesion_count,
        exp_dirs = exp_dirs)
    
    return(adhesion_count_data);
}

gather_general_dynamic_props <- function(results, min.longevity=NA, max.longevity=NA, 
    max.FA.angle=NA,min.FA.angle=NA, max.CHull.percentile = NA, debug=FALSE) {
	points = list()
	for (i in 1:length(results)) {
		res = results[[i]]
        
        if (debug) {
            print(paste("Working on", i))
        }

        filt = ! res$exp_props$split_birth_status & res$exp_props$death_status
        
        if (! is.na(min.longevity)) {
            filt = filt & !is.na(res$exp_props$longevity) & res$exp_props$longevity >= min.longevity;
        }
        if (! is.na(max.longevity)) {
            filt = filt & !is.na(res$exp_props$longevity) & res$exp_props$longevity <= max.longevity;
        }

        if (! is.na(min.FA.angle)) {
            filt = filt & res$exp_props$Mean_FA_recentered_angle > min.FA.angle;
        }
        if (! is.na(max.FA.angle)) {
            filt = filt & res$exp_props$Mean_FA_recentered_angle <= max.FA.angle;
        }

        if (! is.na(max.CHull.percentile)) {
            CHull_dists = res$exp_props$Mean_FA_CHull_dist;
            max_dist = quantile(CHull_dists,max.CHull.percentile)

            filt = filt & CHull_dists <= max_dist;
        }

		points$longevity = c(points$longevity, res$exp_props$longevity[filt])
		points$mean_area = c(points$mean_area, res$exp_props$mean_area[filt])
		points$mean_axial_ratio = c(points$mean_axial_ratio, res$exp_props$mean_axial_ratio[filt])
		points$mean_major_axis = c(points$mean_major_axis, res$exp_props$mean_major_axis[filt])
		points$mean_minor_axis = c(points$mean_minor_axis, res$exp_props$mean_minor_axis[filt])
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

#######################################
# Birth Death Rates
#######################################
count_births_per_point <- function(exp_data) {
    live_adhesions = !is.na(exp_data)
    
    counts = c()
    data_size = dim(live_adhesions)
    for (i in 2:data_size[2]) {
        before_col = live_adhesions[,i-1]
        current_col = live_adhesions[,i]

        birth_rows = which(! before_col & current_col);
        counts = c(counts,length(birth_rows))
    }
    return(counts)
}

count_deaths_per_point <- function(exp_data) {
    live_adhesions = !is.na(exp_data)
    
    counts = c()
    data_size = dim(live_adhesions)
    for (i in 2:data_size[2]) {
        before_col = live_adhesions[,i-1]
        current_col = live_adhesions[,i]

        death_rows = which(before_col & ! current_col);
        counts = c(counts, length(death_rows))
    }
    return(counts)
}

count_live_adhesions_per_point <- function(exp_data) {
    live_adhesions = !is.na(exp_data)

    return(colSums(live_adhesions))
}

matrix_from_list <- function(this_list,lineup='start') {
    max_length = max(unlist(lapply(this_list,length)))

    out_mat = matrix(NA,nrow=length(this_list),ncol=max_length);

    if (lineup == 'start') {
        for (i in 1:length(this_list)) {
            out_mat[i,] = c(this_list[[i]],rep(NA,max_length - length(this_list[[i]])));
        }
    } else if (lineup == 'end') {
        for (i in 1:length(this_list)) {
            out_mat[i,] = c(rep(NA,max_length - length(this_list[[i]])),this_list[[i]]);
        }
    }
    
    return(out_mat)
}

colMedians <- function(this_mat,na.rm=T) {
    medis = c()
    for (i in 1:dim(this_mat)[2]) {
        medis = c(medis,median(this_mat[,i],na.rm=na.rm))
    }
    return(medis)
}

colGeoMeans <- function(this_mat, na.rm=T) {
    geoMeans = c()
    for (i in 1:dim(this_mat)[2]) {
        data = na.omit(this_mat[,i]);
        geoMeans = c(geoMeans,exp(mean(log(data))))
    }
    return(geoMeans)
}

colConfUpper <- function(this_mat) {
    upper = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            upper = c(upper,NA);
        } else {
            #remove any data values that aren't finite, NA, NaN or Inf typically
            data_col = this_mat[,i];
            data_col = data_col[is.finite(data_col)]
            temp = t.test(data_col,conf.level=0.95)
            if (is.nan(temp$conf.int[2])) {
                upper = c(upper,0)
            } else {
                upper = c(upper,temp$conf.int[2])
            }
        }
    }
    return(upper)
}

colConfLower <- function(this_mat) {
    lower = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            lower = c(lower,NA);
        } else {
            #remove any data values that aren't finite, NA, NaN or Inf typically
            data_col = this_mat[,i];
            data_col = data_col[is.finite(data_col)]
            temp = t.test(data_col,conf.level=0.95)
            if (is.nan(temp$conf.int[1])) {
                lower = c(lower,0)
            } else {
                lower = c(lower,temp$conf.int[1])
            }
        }
    }
    return(lower)
}

colSD <- function(this_mat) {
    st_devs = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            st_devs = c(st_devs,NA);
        } else {
            st_devs = c(st_devs,sd(this_mat[,i],na.rm=T))
        }
    }
    return(st_devs)
}

#######################################
# Split Property Measurements
#######################################
count_roll_mean_with_split <- function(exp_data,split.time,exp.norm=T,fold.change=T,
    min.lifetime=NA,dist.measurements = NA, min.dist.ptile = 0.4,use.geo.mean=F) {
    library(zoo)
    
    #######################################################
    # Filtering
    #######################################################
    #check to make sure min.lifetime has been set and that the dimensions of
    #the data are available
    if (! is.na(min.lifetime) & !is.null(dim(exp_data))) {
       lifetimes = rowSums(!is.na(exp_data))
       before_size = dim(exp_data)[1]
       exp_data = exp_data[lifetimes >= min.lifetime,];
       # print(dim(exp_data)[1]/before_size)
    }

    if (! is.na(dist.measurements[1])) {
       before_size = dim(exp_data)[1]
       min.dist = quantile(dist.measurements, c(min.dist.ptile));
       exp_data = exp_data[dist.measurements >= min.dist,]
       # print(dim(exp_data)[1]/before_size)
    }
    
    #######################################################
    # Mean Calc/Splitting/Normalization
    #######################################################
    #dim returns null if the exp_data variable is a single dimension, in that
    #case we want to just use the provided variable as is without taking the
    #column mean
    if (!is.null(dim(exp_data))) {
        if (use.geo.mean) {
            average_val = colGeoMeans(exp_data);
        } else {
            average_val = colMeans(exp_data,na.rm=T);
        }
    } else {
        average_val = exp_data;
    }
    
    average_val = rollmean(average_val,20,na.pad=T);

    average_val_before = average_val[1:split.time]
    average_val_after = average_val[-1:-split.time]

    # average_val_before = rollmean(average_val_before,20,na.pad=T,align='right')
    # average_val_after = rollmean(average_val_after,20,na.pad=T,align='left')
    
    if (exp.norm) {
        mean_before = mean(average_val_before,na.rm=T)
        if (fold.change) {
            average_val_before = average_val_before/mean_before
            average_val_after = average_val_after/mean_before
        } else {
            average_val_before = average_val_before-mean_before
            average_val_after = average_val_after-mean_before
        }
    }

    averages = list()
    averages$before = average_val_before
    averages$after = average_val_after

    return(averages)
}

summarize_split_matrix <- function(exp_data,debug=F) {
    before_temp = matrix_from_list(exp_data$before,lineup='end')
    before_temp = t(tail(t(before_temp),20))
    
    after_temp = matrix_from_list(exp_data$after)
    after_temp = t(head(t(after_temp),70))
    
    results = list()

    results$mat$before = before_temp
    results$mat$after = after_temp

    results$mean$before = colMeans(before_temp,na.rm=T)
    results$mean$after = colMeans(after_temp,na.rm=T)

    results$upper_conf$before = colConfUpper(before_temp)
    results$upper_conf$after = colConfUpper(after_temp)
    
    results$lower_conf$before = colConfLower(before_temp)
    results$lower_conf$after = colConfLower(after_temp)

    results$upper_sd$before = colSD(before_temp) + results$mean$before
    results$upper_sd$after = colSD(after_temp) + results$mean$after
    
    results$lower_sd$before = colSD(before_temp) + results$mean$before
    results$lower_sd$after = colSD(after_temp) + results$mean$after
    return(results);
}

find_mean_split_prop <- function(prop, lineage_props, split.time, min.area = NA, min.data.points=5) {
    prop_size = dim(prop)
    
    prop_before = prop[,1:split.time]
    prop_after = prop[,(split.time+1):prop_size[2]]
    
    #filter out adhesions where we don't see the full lifecycle either before
    #or after drug addition
    before_filter = is.na(prop_before[,1]) & is.na(prop_before[,length(prop_before)])
    after_filter = is.na(prop_after[,1]) & is.na(prop_after[,length(prop_after)])
    
    #filter out adhesions without enough time points
    before_filter = before_filter & rowSums(!is.na(prop_before)) >= min.data.points
    after_filter = after_filter & rowSums(!is.na(prop_after)) >= min.data.points
    
    #apply area filter if requested
    if (!is.na(min.area)) {
        before_filter = before_filter & lineage_props$mean_area >= min.area
        after_filter = after_filter & lineage_props$mean_area >= min.area
    }

    #apply filters
    prop_before = prop_before[before_filter,]
    prop_after = prop_after[after_filter,]

    means = list()
    means$before = rowMeans(prop_before,na.rm=T)
    means$after = rowMeans(prop_after,na.rm=T)
    
    means$area_before = lineage_props$mean_area[before_filter]
    means$area_after = lineage_props$mean_area[after_filter]

    means$longev_before = lineage_props$longevity[before_filter]
    means$longev_after = lineage_props$longevity[after_filter]
    
    return(means)
}

find_cross_split_mean_prop <- function(prop, lineage_props, split.time, area,
    min.area=NA, min.data.points=5) {
    
    #build a filter that only includes adhesions that cross the drug addition
    #time
    all_filter = !is.na(prop[,split.time]) & !is.na(prop[,split.time+1])
    
    prop_size = dim(prop)
    prop_before = prop[,1:split.time]
    prop_after = prop[,(split.time+1):prop_size[2]]
    
    #a minimum length filters to build the before and after addition filters
    before_filter = all_filter & rowSums(!is.na(prop_before)) >= min.data.points
    after_filter = all_filter & rowSums(!is.na(prop_after)) >= min.data.points
    
    #apply area filter if requested
    if (!is.na(min.area)) {
        before_filter = before_filter & lineage_props$mean_area >= min.area
        after_filter = after_filter & lineage_props$mean_area >= min.area
    }
    
    #apply filters
    prop_before = prop_before[before_filter & after_filter,]
    prop_after = prop_after[before_filter & after_filter,]
    
    #split the areas
    area_before = area[,1:split.time]
    area_after = area[,(split.time+1):prop_size[2]]
    area_before = area_before[before_filter & after_filter,]
    area_after = area_after[before_filter & after_filter,]

    means = list()
    means$before = rowMeans(prop_before,na.rm=T)
    means$after = rowMeans(prop_after,na.rm=T)
    means$diff = means$after-means$before
    
    means$area_before = rowMeans(area_before,na.rm=T)
    means$area_after = rowMeans(area_after,na.rm=T)
    means$area_diff = means$area_after-means$area_before

    return(means)
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
	boot_conf = boot.ci(boot_samp,type="bca")
	conf_ints[1,3:4] = boot_conf$perc[4:5]
	conf_ints[1,2] = boot_conf$t0
	bar_lengths[1,1] = boot_conf$t0
	counts[1,1] = length(results_1$a$length)
	if (debug) {
		print('Done 1');
	}
	
	boot_samp_2 = boot(results_2$a$length, function(data,indexes) mean(data[indexes], na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp_2,type="bca")
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
	boot_conf = boot.ci(boot_samp,type="bca")
	conf_ints[2,3:4] = boot_conf$perc[4:5]
	conf_ints[2,2] = boot_conf$t0
	bar_lengths[2,1] = boot_conf$t0
	counts[2,1] = length(results_1$joint$stable_lifetime)
	if (debug) {
		print('Done 3');
	}

	boot_samp_2 = boot(results_2$joint$stable_lifetime, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp_2,type="bca")
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
	boot_conf = boot.ci(boot_samp,type="bca")
	conf_ints[3,3:4] = boot_conf$perc[4:5]
	conf_ints[3,2] = boot_conf$t0
	bar_lengths[3,1] = boot_conf$t0
	counts[3,1] = length(results_1$d$length)
	if (debug) {
		print('Done 5');
	}

	boot_samp_2 = boot(results_2$d$length, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep)
	boot_conf = boot.ci(boot_samp_2,type="bca")
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

gather_barplot_properties <- function(data_sets, bootstrap.rep = 10000) {
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

ranges_overlap <- function(range_1, range_2) {
	if(length(range_1) != 2) browser('range check'); 
	if(length(range_2) != 2) browser('range check');
	stopifnot(length(range_1) == 2)
	stopifnot(length(range_2) == 2)
	
	if (range_1[1] < range_2[1] & range_1[2] > range_2[1]) {
		return(TRUE);
	} else if (range_1[1] > range_2[1] & range_1[1] < range_2[2]) {
		return(TRUE);
	}
	return(FALSE);
}

count_adhesions_per_image <- function(results,time.spacing=1,time.unit=10) {
    data = list();
	for (i in 1:length(results)) {
		res = results[[i]];
        
        num_time_units = (dim(res$exp_data)[[2]]*time.spacing)/time.unit
        num_adhesions = dim(res$exp_data)[[1]]

        data$counts = c(data$counts, num_adhesions/num_time_units);
        data$exp_num = c(data$exp_num, i);
	}

    data$exp_num = as.factor(data$exp_num);
    data = as.data.frame(data);

    return(data);
}

################################################################################
# File Reading/Writing Functions
################################################################################
load_results <- function(dirs,file, debug=TRUE) {
	results = list()
    if (length(dirs) == 0) return()
    if (debug) print(paste('Loading data from',length(dirs),'files.'))

    files_found = 0;
	for (i in 1:length(dirs)) {
		this_file = file.path(dirs[i],file)
		if (file.exists(this_file)) {
            files_found = files_found + 1;
            if (debug) {
                print(paste('Loading File:', this_file))
            }
			var_name = load(file.path(dirs[i],file))
			results[[files_found]] = get(var_name)
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

output_phase_lengths_from_filtered <- function(signif_data,raw_data, dirs = NA) {
    #assembly lengths
    exp_nums = as.numeric(unique(signif_data$assembly$exp_num))
    print('Writing assembly row data to:')
    for (this_exp_num in exp_nums) {
        if (! is.na(dirs)) {
            this_exp_dir = file.path(dirs[this_exp_num],'..');
        } else {
            this_exp_dir = raw_data$intensity[[this_exp_num]]$exp_dir
        }
        print(this_exp_dir)
        this_signif_data = subset(signif_data$assembly, exp_num == this_exp_num, 
            select = c('lin_num','length'))
        write.table(this_signif_data,file=file.path(this_exp_dir,'signif_assembly_rows_lengths.csv'), 
            sep=',', row.names=FALSE, col.names=FALSE)
    }
    
    print('Writing disassembly row data to:')
    exp_nums = as.numeric(unique(signif_data$dis$exp_num))
    for (this_exp_num in exp_nums) {
        if (! is.na(dirs)) {
            this_exp_dir = file.path(dirs[this_exp_num],'..');
        } else {
            this_exp_dir = raw_data$intensity[[this_exp_num]]$exp_dir
        }
        print(this_exp_dir)
        this_signif_data = subset(signif_data$dis, exp_num == this_exp_num, 
            select = c('lin_num','length'))
        write.table(this_signif_data,file=file.path(this_exp_dir,'signif_disassembly_rows_lengths.csv'), 
            sep=',', row.names=FALSE, col.names=FALSE)
    }
}

###############################################################################
#Spacial Functions
###############################################################################
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

bin_corr_data <- function(corr_results, bin_size = NA, bootstrap.rep = 5000, 
    bin_max = NA, pixel_size=0.215051) {
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
			boot_conf = boot.ci(boot_samp,type="bca")
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
