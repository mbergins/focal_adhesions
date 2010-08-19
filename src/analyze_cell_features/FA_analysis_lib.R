################################################################################
# FA_analysis_lib.R: various functions used to find and plot the linear regions
#  and associated data from the focal adhesion identification/analysis programs
################################################################################

########################################
#Data fitting functions
########################################
gather_bilinear_models_from_dirs <- function (dirs, min_length = 10,
        data_file='Average_adhesion_signal.csv', col_lims = NA, 
        normed = TRUE, log.trans = TRUE, boot.samp = NA, results.file = NA,
        debug = FALSE) {

    results = list();

    for (i in 1:length(dirs)) {
        if (is.na(dirs[[i]])) {
            next;
        }
        if (! file.exists(file.path(dirs[[i]],data_file))) {
            next;
        }
        
        if (debug) {
            print(dirs[[i]])
        }

        exp_data <- as.matrix(read.table(file.path(dirs[[i]],data_file),header = FALSE, sep = ','));

        exp_props <- read.table(file.path(dirs[[i]],'../single_lin.csv'), header = TRUE, sep = ',');

        #process the col_lim parameter passed in if values were passed in
        this_col_lim = NA;
        if (! is.na(as.matrix(col_lims)[1,1])) {
            if(dim(as.matrix(col_lims))[[2]] == 1) {
                this_col_lim = c(col_lims[i],dim(exp_data)[[2]])
            } else {
                this_col_lim = col_lims[i,]
            }
        }

        results[[i]] <- gather_bilinear_models(exp_data, exp_props, 
                min_length = min_length, col_lims = this_col_lim, 
                normed = normed, log.trans = log.trans, 
                boot.samp = boot.samp, debug=debug);

        regex_range = regexpr("time_series_[[:digit:]].*?/",dirs[[1]], perl=TRUE);
        if (regex_range[1] == -1) {
            results[[i]]$exp_dir = dirs[[i]];
        } else {
            results[[i]]$exp_dir = substr(dirs[[i]], regex_range[1], regex_range[1] + attr(regex_range,'match.length') - 2);
        }
        results[[i]]$parameters$normed = normed;
        results[[i]]$parameters$min_length = min_length;
        results[[i]]$parameters$log.trans = log.trans;

        if (! is.na(results.file)) {
            this_file = file.path(dirs[[i]],results.file);
            if (! file.exists(dirname(this_file))) {
                dir.create(dirname(this_file),recursive=TRUE);
            }
            this_result = results[[i]];
            save(this_result,file = this_file);
        }
    }

    return(results);
}

gather_bilinear_models <- function(data_set, props, 
    min_length = 10, col_lims = NA, normed = TRUE, 
    log.trans = TRUE, boot.samp = NA, debug = FALSE) {
        
    if (is.numeric(col_lims) && length(col_lims) == 2) {
        data_set = data_set[,col_lims[1]:col_lims[2]];
    }

    rows <- dim(data_set)[[1]];
    cols <- dim(data_set)[[2]];

    results <- list();
    results$exp_props = props;

    results$exp_data = data_set;

    results$stable_data_set = list();
    for (i in 1:rows) {
        if (i > 100 & debug) {
            #next();
        }

        if ((i %% 100 == 0 | i < -300) & debug) {
            print(sprintf('%d %d %.02f',i,rows,i/rows))
        }

        this_data_set = as.vector(data_set[i,]);
        these_exp_props = results$exp_props[i,];

        temp_results = find_optimum_bilinear_fit(this_data_set, these_exp_props, normed = normed, 
                min_length = min_length, log.trans = log.trans);

        if (is.na(charmatch("assembly", names(results)))) {
            results$assembly = temp_results$assembly;
        } else {
            stopifnot(length(temp_results$assembly) == dim(results$assembly)[[2]]);
            results$assembly = rbind(results$assembly, temp_results$assembly);
        }
        if (dim(results$assembly)[[1]] != i) {
            print(paste("Problem with row count in row number:", i));
            stop();
        }

        if (is.na(charmatch("disassembly", names(results)))) {
            results$disassembly = temp_results$disassembly;
        } else {
            stopifnot(length(temp_results$disassembly) == dim(results$disassembly)[[2]]);
            results$disassembly = rbind(results$disassembly, temp_results$disassembly);
        }
        stopifnot(dim(results$disassembly)[[1]] == i);

        #If either of the length values are NA or the assembly and disassembly phases take up the entire data set, 
        if (! is.na(results$disassembly$length[i]) && ! is.na(results$assembly$length[i]) &&
            ! is.na(results$exp_props$longevity[i]) &&
            results$exp_props$longevity[i] > (results$assembly$length[i] + results$disassembly$length[i])) {

            numeric_data_set = na.omit(this_data_set);

            results$stable_data_set[[i]] = numeric_data_set[(results$assembly$length[i]+1):(length(numeric_data_set) - results$disassembly$length[i])];
            results$stable_lifetime[i] = length(results$stable_data_set[[i]]);
            results$stable_mean[i] = mean(results$stable_data_set[[i]]);
            results$stable_variance[i] = var(results$stable_data_set[[i]]);

            stopifnot(length(numeric_data_set) == results$stable_lifetime[i] + results$assembly$length[i] + results$disassembly$length[i]);
        }
    }

    if (debug) {
        print('Done with gathering bilinear models-Starting to Pad Results');
    }
    results <- pad_results_to_row_length(results, rows, debug=debug);
    if (debug) {
        print('Done with Padding Results');
    }

    if (is.numeric(boot.samp)) {
        results$sim_results <- gather_linear_regions.boot(results, min_length = min_length, 
                col_lims = col_lims, normed = normed, log.trans = log.trans, boot.samp = boot.samp);
    }
        print(names(results));
    results
}

pad_results_to_row_length <- function(results, desired_length, debug=FALSE) {
	print(dim(results$assembly)[[1]])
	print(desired_length)
	stopifnot(dim(results$assembly)[[1]] <= desired_length)
	while (dim(results$assembly)[[1]] != desired_length) {
		results$assembly = rbind(results$assembly, rep(NA, length(results$assembly[1,])));
	}
	stopifnot(desired_length == dim(results$assembly)[[1]])
	if (debug) {
		print('Done with assembly padding')
	}
	
	stopifnot(dim(results$disassembly)[[1]] <= desired_length)
	while (dim(results$disassembly)[[1]] != desired_length) {
		results$disassembly = rbind(results$disassembly, rep(NA, length(results$disassembly[1,])));
	}
	if (debug) {
		print('Done with disassembly padding')
	}

	stopifnot(length(results$stable_lifetime) == length(results$stable_data_set))
	stopifnot(length(results$stable_lifetime) == length(results$stable_variance))
	stopifnot(length(results$stable_lifetime) == length(results$stable_mean))
		
	stopifnot(length(results$stable_lifetime) <= desired_length)
	while (length(results$stable_lifetime) < desired_length) {
		results$stable_lifetime = c(results$stable_lifetime, NA);
		results$stable_data_set = c(results$stable_data_set, NA);
		results$stable_variance = c(results$stable_variance, NA);
		results$stable_mean = c(results$stable_mean, NA);
	}
	if (debug) {
		print('Done with stability padding')
	}

	results
}

find_optimum_bilinear_fit <- function(initial_data_set, exp_props, 
	normed = TRUE, min_length = 10, log.trans = TRUE) {

	results = list()
	resid = list(assembly = list(), disassembly = list())
	
	results$filt_init = initial_data_set[! is.nan(initial_data_set)]
	if (length(results$filt_init) == 0) {
		stop('Length of data set is equal to zero, check to make sure each line has data in it')
	}
	
	this_data_set = data.frame(y = results$filt_init, x = 1:length(results$filt_init))
	#When building models with the overall cell background corrected intensity data, there 
	#is occasionally an entry in a lineage that is below zero. Taking the log of 
	#a negative number causes problems, so I'll search for negative entries in the 
	#lineage and replace them with the smallest number in the data set that is greater 
	#than zero. If none of the data entries are greater than zero, set the negative 
	#entires to a small number.
	if (length(this_data_set$y[this_data_set$y > 0]) == 0) {
		min_entry = 0.00001;
	} else {
		min_entry = min(this_data_set$y[this_data_set$y > 0])
	}
	this_data_set$y[this_data_set$y <= 0] = min_entry;
	
	#Search the beginning of the sequence for a linear fit
	assembly_slope_calculated = FALSE;
    # we don't want to find an assembly slope if the first entry in the raw
    # data set isn't a "NA", as that indicates that we saw the birth of the
    # adhesion during the experiment, we are also not interested in adhesions
    # that are the result of a split birth event, finally, we exclude those
    # adhesions that aren't long enough to accomodate the minimum assembly and
    # disassembly phases
	if (is.nan(initial_data_set[1]) & ! exp_props$split_birth_status & length(results$filt_init) >= (min_length * 2)) {
		assembly_slope_calculated = TRUE;
		for (j in min_length:dim(this_data_set)[[1]]) {
			assembly_subset = this_data_set[1:j,]
			stopifnot(length(assembly_subset$y) >= min_length)
			stopifnot(length(assembly_subset$y) == j)
			if (normed) {
				assembly_subset$y = assembly_subset$y/assembly_subset$y[1]
			}
			if (log.trans) {
				assembly_subset$y = log(assembly_subset$y)
			}

			model <- lm(y ~ x, data = assembly_subset)
			summary <- summary(model);
			
			results$assembly$R_sq[j] = summary$adj.r.squared
			#dealing with a degerate case, where lm produce NaN for the R squared 
			#value when the data set is a flat line, see:
			#	>data <- data.frame(x = c(1,2,3), y = c(1,1,1))
			#	>summary(lm(y ~ x, data=data))
			if (is.nan(results$assembly$R_sq[j])) {
				results$assembly$R_sq[j] = 1
			}
			
			results$assembly$p_val[j] = summary$coefficients[2,4]
			results$assembly$length[j] = dim(assembly_subset)[[1]]
			results$assembly$inter[j] = coef(model)[[1]]
			results$assembly$slope[j] = coef(model)[[2]]
			results$assembly$first_point[j] = assembly_subset$y[1]
			results$assembly$fold_change[j] = max(assembly_subset$y)/min(assembly_subset$y)
			
			resid$assembly[[j]] = as.numeric(resid(model))
		}
	} else {
		results$assembly$R_sq[1] = 0

		results$assembly$p_val[1] = NA		
		results$assembly$length[1] = NA
		results$assembly$inter[1] = NA
		results$assembly$slope[1] = NA
		results$assembly$first_point[1] = NA
		results$assembly$fold_change[1] = NA
	
		resid$assembly[[1]] = NA
	}
	
	#Search the end of the sequence for a linear fit
	disassembly_slope_calculated = FALSE;
    # we don't want to find a disassembly slope if the last entry in the raw
    # data set isn't a "NA", as that indicates that we didn't see the death of
    # the adhesion during the experiment, we are also not interested in
    # adhesions that merge another adhesion, finally, we exclude those
    # adhesions that aren't long enough to accomodate the minimum assembly and
    # disassembly phases
	if (is.nan(initial_data_set[length(initial_data_set)]) & exp_props$death_status & length(results$filt_init) >= (min_length * 2)) {
		disassembly_slope_calculated = TRUE;
		for (j in min_length:dim(this_data_set)[[1]]) {
			disassembly_subset = this_data_set[(dim(this_data_set)[[1]]-j+1):dim(this_data_set)[[1]],]
			
			if (normed) {
				disassembly_subset$y = disassembly_subset$y[1]/disassembly_subset$y
			}
			if (log.trans) {
				disassembly_subset$y = log(disassembly_subset$y)
			}
	
			model <- lm(y ~ x, data = disassembly_subset)
			summary <- summary(model);
			
			results$disassembly$R_sq[j] = summary$adj.r.squared
			#dealing with a degerate case, where lm produces NaN for the R squared 
			#value when the data set is a flat line, see:
			#	>data <- data.frame(x = c(1,2,3), y = c(1,1,1))
			#	>summary(lm(y ~ x, data=data))
			if (is.nan(results$disassembly$R_sq[j])) {
				results$disassembly$R_sq[j] = 1
			}
			
			results$disassembly$p_val[j] = summary$coefficients[2,4]
			results$disassembly$length[j] = dim(disassembly_subset)[[1]]
			results$disassembly$inter[j] = coef(model)[[1]]
			results$disassembly$slope[j] = coef(model)[[2]]
			results$disassembly$first_point[j] = disassembly_subset$y[1]
			results$disassembly$fold_change[j] = max(disassembly_subset$y)/min(disassembly_subset$y)
			resid$disassembly[[j]] = as.numeric(resid(model))
		}
	} else {
		results$disassembly$R_sq[1] = 0
		
		results$disassembly$p_val[1] = NA	
		results$disassembly$length[1] = NA
		results$disassembly$inter[1] = NA
		results$disassembly$slope[1] = NA
		results$disassembly$first_point[1] = NA
		results$disassembly$fold_change[1] = NA
		
		resid$disassembly[[1]] = NA	
	}
	
	best_indexes = c()
	if (! any(assembly_slope_calculated, disassembly_slope_calculated)) {
		results$assembly$R_sq[1] = NA
		results$disassembly$R_sq[1] = NA
		best_indexes = c(1,1);
	} else {
		best_indexes = find_best_length_combination(results, min_length = min_length)
	
		#With the R squared matrix calculated reset the r_sq componenets to NA, if needed since 
		#there were no fits calculated for them
		if (! assembly_slope_calculated) {
			results$assembly$R_sq[1] = NA
		}	
		if (! disassembly_slope_calculated) {
			results$disassembly$R_sq[1] = NA
		}
	}
	
	best_results = list()
				
	best_results$assembly = as.data.frame(results$assembly)[best_indexes[1],]
	best_results$assembly$residual = resid$assembly[best_indexes[1]]
	best_results$disassembly = as.data.frame(results$disassembly)[best_indexes[2],]
	best_results$disassembly$residual = resid$disassembly[best_indexes[2]]

	best_results
}

find_best_length_combination <- function(results, min_length = 10) {
	#Build an array with the sums of the collected R square values
	R_sq_sums = array(NA, c(length(results$assembly$R_sq),length(results$disassembly$R_sq)));
	for (i in 1:length(results$assembly$R_sq)) {
		if (is.na(results$assembly$R_sq[i]) | is.nan(results$assembly$R_sq[i])) {
			next
		}
		for (j in 1:length(results$disassembly$R_sq)) {
			if (is.na(results$disassembly$R_sq[j]) | is.nan(results$disassembly$R_sq[j])) {
				next
			}
			if ((j+i) > length(results$filt_init)) {
				next
			}
			
			R_sq_sums[i,j] = results$assembly$R_sq[i]+results$disassembly$R_sq[j]
		}
	}
		
	#locate positions of the highest R squared value
	max_R_sq = max(R_sq_sums[! is.na(R_sq_sums)])
	if (max_R_sq == -Inf) {
		print(R_sq_sums)
	}
	
	highest_square_priority = -Inf
	
	for (i in 1:dim(R_sq_sums)[[1]]) {
		for (j in 1:dim(R_sq_sums)[[2]]) {
			if (is.na(R_sq_sums[i,j]) | is.nan(R_sq_sums[i,j])) {
				next
			}
			
			i_priority = i - min_length + 1
			if (i_priority <= 0) {
				i_priority = 0
			}
			j_priority = j - min_length + 1
			if (j_priority <= 0) {
				j_priority = 0
			}
			this_priority = i_priority + j_priority
			
			if (  max_R_sq == R_sq_sums[i,j]
				& highest_square_priority < this_priority) {
				
				highest_square_priority = this_priority
				best_indexes = c(i,j)
			}
		}
	}
	
	best_indexes	
}

gather_linear_regions.boot <- function(results, 
	min_length = 10, col_lims = NaN, normed = 1, 
	log.trans = TRUE, boot.samp = NA) {

	sim_results <- list()
	
        #collect the entire set of adhesion signal values and lengths,
        #excluding adhesions which don't live long enough
	
        all_ad_sig = c()
	all_length = c()
	for (i in 1:dim(results$exp_data)[[2]]) {
		temp = as.numeric(results$exp_data[i,])
		temp = temp[! is.nan(temp)]
		if (length(temp) >= min_length * 2) {
			all_length = c(all_length,length(temp))
			all_ad_sig = c(all_ad_sig,temp)	
		}
	}
	
	#produce fake adhesion signal and props variables by sampling from the values collected above
	sim_ad_sig <- array(NaN, dim = c(boot.samp, max(all_length) + 2))
	sim_props = list()
	for (i in 1:boot.samp) {
		this_length = sample(all_length,1)
		data = sample(all_ad_sig, this_length, replace = TRUE)
		sim_ad_sig[i,] <- c(NaN, data, array(NaN, dim = c(max(all_length) - length(data) + 1)))
		sim_props$death_status[i] = 1
	}
	
	sim_results <- gather_bilinear_models(sim_ad_sig, sim_props, 
				       min_length = min_length, normed = normed, 
                                       log.trans = log.trans, save.exp_data = FALSE)
        sim_results
}	

gather_correlations_from_dirs <- function (dirs, results, data_file='Area.csv',
	result.normed = TRUE, exp.normed = TRUE, result.log.trans = TRUE, 
	exp.log.trans = TRUE, results.file = NA, save.exp_data = TRUE) {
	
	corr_results = list()
	
	for (k in 1:length(dirs)) {
		if (is.na(dirs[[k]])) {
			next
		}
		
		print(dirs[[k]])
		
		exp_data <- as.matrix(read.table(file.path(dirs[[k]],data_file),header = FALSE, sep  = ','));
		
		corr_results[[k]] <- gather_correlations(results[[k]], exp_data, 
			result.normed = result.normed, exp.normed = exp.normed, 
			result.log.trans = result.log.trans, exp.log.trans = exp.log.trans, 
			results.file = results.file, save.exp_data = save.exp_data);
			
                if (! is.na(results.file)) {
                    this_result = corr_results[[k]]
                        save(this_result,file = file.path(dirs[[k]],results.file))
                }
	}
	
	corr_results
}

gather_correlations <- function(result, exp_data, result.normed = TRUE, 
	exp.normed = FALSE, result.log.trans = TRUE, exp.log.trans = FALSE, 
	results.file = NA, save.exp_data = TRUE) {

	corr_result = list()
	
	if (save.exp_data) {
		corr_result$exp_data = exp_data
	}
	
	count = 0
	for (i in 1:length(result$assembly$R_sq)) {
		data_1 = as.numeric(result$exp_data[i,])
		data_1 = data_1[! is.nan(data_1)]
		data_2 = as.numeric(exp_data[i,])
		data_2 = data_2[! is.nan(data_2)]		
		
		corr_result$assembly[i] = NA
		corr_result$disassembly[i] = NA

		if (! is.na(result$assembly$R_sq[i])) {
			this_data_1 = data_1[1:result$assembly$length[i]]
			this_data_2 = data_2[1:result$assembly$length[i]]
			
			if (result.normed) {
				this_data_1 = this_data_1/this_data_1[1]
			}
			if (exp.normed) {
				this_data_2 = this_data_2/this_data_2[1]
			}
			
			if (result.log.trans) {
				this_data_1 = log(this_data_1)
			}
			if (exp.log.trans) {
				this_data_2 = log(this_data_2)
			}
			
			if (sd(this_data_1) != 0 & sd(this_data_2) != 0) {
				corr_result$assembly[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$estimate)
				corr_result$conf$assembly_lower[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[1]
				corr_result$conf$assembly_upper[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[2]
			}
		}
		if (! is.na(result$disassembly$R_sq[i])) {
			this_data_1 = data_1[(length(data_1) - result$disassembly$length[i]):length(data_1)]
			this_data_2 = data_2[(length(data_2) - result$disassembly$length[i]):length(data_2)]

			if (result.normed) {
				this_data_1 = this_data_1[1]/this_data_1
			}
			if (exp.normed) {
				this_data_2 = this_data_2[1]/this_data_2
			}
			
			if (result.log.trans) {
				this_data_1 = log(this_data_1)
			}
			if (exp.log.trans) {
				this_data_2 = log(this_data_2)
			}
				
			if (sd(this_data_1) != 0 & sd(this_data_2) != 0) {
				corr_result$disassembly[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$estimate)
				corr_result$conf$disassembly_lower[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[1]
				corr_result$conf$disassembly_upper[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[2]
			}
		}
	}
	corr_result
}

########################################
#Plotting Functions
########################################
plot_ad_seq <- function (results,index,type='assembly',log.trans = TRUE, ...) {
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

boxplot_with_points <- function(data,
    colors=c('darkgreen','red','yellow','blue','pink','cyan','gray','orange','brown','purple'),
    notch=F, names, range=1.5, inc.n.counts = TRUE, inc.points = TRUE, pch=20, na.omit = TRUE,
    point_cex=0.5, return_output = FALSE, with.median.props = TRUE, 
    median.props.pos = c(0.5,0.9), median.props.color = 'blue', ...) {
	
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
                   temp_data,col=colors[[i]], pch=pch, cex=point_cex, )
		}
	}
    if (with.median.props) {
        p_vals = determine_median_p_value(data[[1]], data[[2]]);
        median_ratio_low = sprintf('%.0f', 100*p_vals$ratio_conf[1]);
        median_ratio_high = sprintf('%.0f', 100*p_vals$ratio_conf[2]);
        
        median_ratio_average = (p_vals$median_vals[2] - p_vals$median_vals[1])/p_vals$median_vals[1];
        median_ratio_average = sprintf('%.0f', 100*median_ratio_average);

        x_pos = (plot_dims[2] - plot_dims[1])*median.props.pos[1] + plot_dims[1]
        y_pos = (plot_dims[4] - plot_dims[3])*median.props.pos[2] + plot_dims[3]

        text(x_pos,y_pos,
             paste('p<',p_vals$p_val, sep=''), 
             col=median.props.color);
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
filter_results <- function(results, min_R_sq=0.9, max_p_val = 0.05, debug = FALSE,
        pos_slope = TRUE, primary_filter_results = NA, cell_intensities = NA) {

    stopifnot(is.na(primary_filter_results) | length(results) == length(primary_filter_results));

    points = list();
    for (i in 1:length(results)) {
        if (debug) {
            print(paste("Working on Experiment #:",i));
        }
        res = results[[i]];

        if (all(is.na(primary_filter_results))) {
            for_filter = res
        } else {
            for_filter = primary_filter_results[[i]]
        }

        filter_sets = produce_rate_filters(for_filter, min_R_sq = min_R_sq, max_p_val = max_p_val, pos_slope = pos_slope)
        
        assembly_filt = filter_sets$assembly;
        disassembly_filt = filter_sets$disassembly;
        joint_filt = filter_sets$joint;

        points$assembly$slope = c(points$assembly$slope, res$assembly$slope[assembly_filt]);
        points$assembly$R_sq = c(points$assembly$R_sq, res$assembly$R_sq[assembly_filt]);
        points$assembly$p_val = c(points$assembly$p_val, res$assembly$p_val[assembly_filt]);
        points$assembly$length = c(points$assembly$length, res$assembly$length[assembly_filt]);
        
        points$assembly$longevity = c(points$assembly$longevity, res$exp_props$longevity[assembly_filt]);
        points$assembly$start_x = c(points$assembly$start_x, res$exp_props$start_x[assembly_filt]);
        points$assembly$start_y = c(points$assembly$start_y, res$exp_props$start_y[assembly_filt]);
        points$assembly$mean_area = c(points$assembly$mean_area, res$exp_props$mean_area[assembly_filt]);
 
        points$assembly$lin_num = c(points$assembly$lin_num, which(assembly_filt));
        points$assembly$exp_dir = c(points$assembly$exp_dir, rep(res$exp_dir, length(which(assembly_filt))));
        points$assembly$exp_num = c(points$assembly$exp_num, rep(i,length(which(assembly_filt))));
        for (j in which(assembly_filt)) {
            numeric_indexes = which(! is.na(res$exp_data[j,]));
            numeric_indexes = numeric_indexes[1];
            points$assembly$start = c(points$assembly$start, mean(res$exp_data[j,numeric_indexes]));
            
            numeric_indexes = which(! is.na(res$exp_data[j,]));
            numeric_indexes = numeric_indexes[(res$assembly$length[j]) : res$assembly$length[j]];
            points$assembly$end = c(points$assembly$end, mean(res$exp_data[j,numeric_indexes]));
        }
        if (any(names(res$exp_props) == 'starting_edge_dist')) {
            points$assembly$edge_dist = c(points$assembly$edge_dist, res$exp_props$starting_edge_dist[assembly_filt])
            points$control$birth_pos = c(points$control$birth_pos,  mean(res$exp_props$starting_edge_dist[assembly_filt], na.rm=T));
        }
        for (j in names(points$assembly)) {
            if (length(points$assembly[[1]]) != length(points$assembly[[j]])) {
                print(paste("Length of property '",j, "' wrong in assembly filtering", sep=''))
                print(paste("Lengths are", length(points$assembly[[1]]), " ", length(points$assembly[[j]])))
                browser()
                stop()
            }
        }

        points$disassembly$slope = c(points$disassembly$slope, res$disassembly$slope[disassembly_filt]);
        points$disassembly$R_sq = c(points$disassembly$R_sq, res$disassembly$R_sq[disassembly_filt]);
        points$disassembly$p_val = c(points$disassembly$p_val, res$disassembly$p_val[disassembly_filt]);
        points$disassembly$length = c(points$disassembly$length, res$disassembly$length[disassembly_filt]);

        points$disassembly$longevity = c(points$disassembly$longevity, res$exp_props$longevity[disassembly_filt]);
        points$disassembly$start_x = c(points$disassembly$start_x, res$exp_props$start_x[disassembly_filt]);
        points$disassembly$start_y = c(points$disassembly$start_y, res$exp_props$start_y[disassembly_filt]);
        points$disassembly$mean_area = c(points$disassembly$mean_area, res$exp_props$mean_area[disassembly_filt]);
        
        points$disassembly$lin_num = c(points$disassembly$lin_num, which(disassembly_filt));
        points$disassembly$exp_dir = c(points$disassembly$exp_dir, rep(res$exp_dir,length(which(disassembly_filt))));
        points$disassembly$exp_num = c(points$disassembly$exp_num, rep(i,length(which(disassembly_filt))));
        for (j in which(disassembly_filt)) {
            numeric_indexes = which(! is.na(res$exp_data[j,]));
            numeric_indexes = numeric_indexes[(length(numeric_indexes)-res$disassembly$length[j]):length(numeric_indexes)];
            numeric_indexes = numeric_indexes[1]
            points$disassembly$start = c(points$disassembly$start, mean(res$exp_data[j,numeric_indexes]));
            
            numeric_indexes = which(! is.na(res$exp_data[j,]));
            numeric_indexes = numeric_indexes[(length(numeric_indexes)-res$disassembly$length[j]):length(numeric_indexes)];
            numeric_indexes = numeric_indexes[(length(numeric_indexes)):length(numeric_indexes)];
            points$disassembly$end = c(points$disassembly$end, mean(res$exp_data[j,numeric_indexes]));
        }
        if (any(names(res$exp_props) == 'ending_edge_dist')) {
            points$disassembly$edge_dist = c(points$disassembly$edge_dist, res$exp_props$ending_edge_dist[disassembly_filt])
            points$control$death_pos = c(points$control$death_pos,  mean(res$exp_props$ending_edge_dist[disassembly_filt], na.rm=T));
        }
        for (j in names(points$disassembly)) {
            if (length(points$disassembly[[1]]) != length(points$disassembly[[j]])) {
                print(paste("Length of property '",j, "' wrong in disassembly filtering", sep=''))
                print(paste("Lengths are", length(points$disassembly[[1]]), " ", length(points$disassembly[[j]])))
                stop()
            }
        }

        points$joint$birth_dist = c(points$joint$birth_dist, res$exp_props$starting_edge_dist[joint_filt]);
        points$joint$death_dist = c(points$joint$death_dist, res$exp_props$ending_edge_dist[joint_filt]);

        points$joint$lin_num = c(points$joint$lin_num, which(joint_filt));
        points$joint$exp_num = c(points$joint$exp_num, rep(i,length(which(joint_filt))));

        points$joint$assembly = c(points$joint$assembly, res$assembly$slope[joint_filt]);
        points$joint$assembly_R = c(points$joint$assembly_R, res$assembly$R_sq[joint_filt]);
        points$joint$disassembly = c(points$joint$disassembly, res$disassembly$slope[joint_filt]);
        points$joint$disassembly_R = c(points$joint$disassembly_R, res$disassembly$R_sq[joint_filt]);

        points$joint$assembly_length = c(points$joint$assembly_length, res$assembly$length[joint_filt]);
        points$joint$disassembly_length = c(points$joint$disassembly_length, res$dis$length[joint_filt]);
        points$joint$stable_lifetime = c(points$joint$stable_lifetime, res$stable_lifetime[joint_filt]);
        points$joint$total_lifetime = c(points$joint$total_lifetime, res$exp_props$longevity[joint_filt]);
        points$joint$stable_mean = c(points$joint$stable_mean, res$stable_mean[joint_filt]);
        points$joint$stable_variance = c(points$joint$stable_variance, res$stable_variance[joint_filt]);
        for (j in names(points$joint)) {
            if (length(points$joint[[1]]) != length(points$joint[[j]])) {
                print(paste("Length of property '",j, "' wrong in joint filtering", sep=''))
                print(paste("Lengths are", length(points$joint[[1]]), " ", length(points$joint[[j]])))
                stop()
            }
        }
        
        if (! all(is.na(cell_intensities))) {
            points$control$assembly_slope = c(points$control$assembly_slope,  mean(res$assembly$slope[assembly_filt]));
            points$control$disassembly_slope = c(points$control$disassembly_slope,  mean(res$disassembly$slope[disassembly_filt]));
            points$control$mean_cell_int = c(points$control$mean_cell_int, mean(cell_intensities[[i]]));
        }
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

produce_rate_filters <- function(raw_data, min_R_sq=0.9, max_p_val=0.05, pos_slope=TRUE) {
    filter_sets = list()

    # First we cycle through all the variables used to filter the assembly and
    # disassembly phases, collecting logic table describing which adhesions
    # fulfill which criteria
    vars_to_filter = c("assembly", "disassembly")
    for (var in vars_to_filter) {
        filter_sets[[var]]$good_R_sq = ! is.na(raw_data[[var]]$R_sq) & raw_data[[var]]$R_sq >= min_R_sq
        filter_sets[[var]]$low_p_val = ! is.na(raw_data[[var]]$p_val) & raw_data[[var]]$p_val < max_p_val
        filter_sets[[var]]$pos_slope = ! is.na(raw_data[[var]]$slope) & raw_data[[var]]$slope > 0; 
    }
    
    # There are two variables that are unique to assembly and disassembly
    # (split birth events and deaths due to merges), we will keep all those in
    # the logic table "extra" and the in named variables so we can access the
    # named variables separated if needed
    filter_sets$assembly$not_split_birth = ! raw_data$exp_props$split_birth_status
    filter_sets$assembly$extra = filter_sets$assembly$not_split_birth
    
    filter_sets$disassembly$death_status = raw_data$exp_props$death_status
    filter_sets$disassembly$extra = filter_sets$disassembly$death_status
    
    # Now all the filters are cascaded to produce the final filter set, then
    # the pos_slope variable is checked and applied if requested
    final_filters = list()
    for (var in vars_to_filter) {
        final_filters[[var]] = (filter_sets[[var]]$good_R_sq & filter_sets[[var]]$low_p_val & filter_sets[[var]]$extra);
    }

    if (pos_slope) {
        for (var in vars_to_filter) {
            final_filters[[var]] = (final_filters[[var]] & filter_sets[[var]]$pos_slope);
        }
    }

    final_filters$joint = final_filters$assembly & final_filters$disassembly
    
    for (filter_type in names(final_filters)) {
        final_filters[[filter_type]] = as.logical(final_filters[[filter_type]])
    }

    final_filters$filter_sets = filter_sets;
    
    return(final_filters)
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

determine_mean_conf_int <- function(data, bootstrap.rep = 10000) {
	require(boot);
	boot_samp = boot(data, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep);
		
    conf_int = boot.ci(boot_samp, type="perc", conf=0.95)$percent[4:5]

    return(conf_int)
}

determine_mean_p_value <- function(data_1,data_2, bootstrap.rep = 10000) {
	require(boot);
	boot_samp_1 = boot(data_1, function(data_1,indexes) mean(data_1[indexes],na.rm=T), bootstrap.rep);
	boot_samp_2 = boot(data_2, function(data_1,indexes) mean(data_1[indexes],na.rm=T), bootstrap.rep);
		
    conf_int_one = boot.ci(boot_samp_1, type="perc", conf=0.95)
    conf_int_two = boot.ci(boot_samp_2, type="perc", conf=0.95)

    results = list();

    results$p_val = find_p_val_from_bootstrap(boot_samp_1, boot_samp_2);
    results$mean_vals = c(boot_samp_1$t0, boot_samp_2$t0);
    results$conf_ints = rbind(conf_int_one$percent[4:5], conf_int_two$percent[4:5])

    return(results);
}

determine_median_p_value <- function(data_1,data_2, bootstrap.rep = 10000) {
	require(boot);
    
    data_package = list(one = data_1, two = data_2);
    boot_ratio = boot(data_package, function(values, indexes) ratio_samp(values$one, values$two), bootstrap.rep);
    boot_ratio_conf = boot.ci(boot_ratio,type="perc", conf=0.99)

    results = list()
	
    results$ratio_conf = boot_ratio_conf$perc[4:5];
	
    boot_samp_1 = boot(data_1, function(values,indexes) median(values[indexes],na.rm=T), bootstrap.rep);
	boot_samp_2 = boot(data_2, function(values,indexes) median(values[indexes],na.rm=T), bootstrap.rep);
	results$p_val = find_p_val_from_bootstrap(boot_samp_1, boot_samp_2);
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
	p_vals_to_test = c(seq(0.9,0.11,by=-0.01),0.1,0.05,1E-2,1E-3,1E-4,1E-5)) {
	
	stopifnot(class(boot_one) == "boot")
	stopifnot(class(boot_two) == "boot")
	
	overlap_region = c()
	for (this_val in p_vals_to_test) {
		if (any(overlap_region)) {
			next;
		}
		conf_int_one = boot.ci(boot_one, type="perc", conf=(1-this_val))
		conf_int_two = boot.ci(boot_two, type="perc", conf=(1-this_val))
				
		if (ranges_overlap(conf_int_one$perc[4:5], conf_int_two$perc[4:5])) {
			overlap_region = c(overlap_region, TRUE)
		} else {
			overlap_region = c(overlap_region, FALSE)
		}
	}
	
	overlap_indexes = which(overlap_region)
	if (length(overlap_indexes) == 0) {
		return(p_vals_to_test[length(p_vals_to_test)])
	} else if (overlap_indexes[1] == 1) {
		return(overlap_indexes[1]);
	} else {
		return(p_vals_to_test[overlap_indexes[1] - 1])
	}	
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
		points$ending_edge = c(points$ending_edge, res$exp_props$ending_edge[filt])
		points$starting_edge = c(points$starting_edge, res$exp_props$starting_edge[filt])
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
	for (i in 1:length(dirs)) {
		this_file = file.path(dirs[i],file)
		if (file.exists(this_file)) {
            if (debug) {
                print(paste('Loading File:', this_file))
            }
			load(file.path(dirs[i],file))
			results[[i]] = this_result
		}
	}
	results
}

load_results_2 <- function(dirs,file) {
	temp_results = list()
	for (i in 1:length(dirs)) {
		this_file = file.path(dirs[i],file)
		if (file.exists(this_file)) {
			load(file.path(dirs[i],file))
			temp_results[[i]] = as.data.frame(results)
		}
	}
	if (length(temp_results) == 1) {
		temp_results = results[[1]]
	}
	temp_results
}

load_results_data_frame <- function(dirs,file,variable_name, debug=FALSE) {
    #simple function to read in all the R data files in each of the specified
    #dirs, searching for the provided variable name after every load, note that
    #all the data is concatenated into a single data frame
	
    results = data.frame()
	for (i in 1:length(dirs)) {
		this_file = file.path(dirs[i],file)
        if (debug) {
            print(paste("Loading:", file.path(dirs[i],file)))
        }
		if (file.exists(this_file)) {
			load(file.path(dirs[i],file))
			results = rbind(results,get(variable_name))
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

write_assembly_disassembly_periods <- function(result, dir) {
	if (! file.exists(dir)) {
		dir.create(dir,recursive=TRUE)
	}
	
	row_nums = which(! is.na(result$assembly$length))
	if (! is.null(row_nums)) {
		rows_and_length = cbind(row_nums, result$assembly$length[row_nums]);
		write.table(rows_and_length,file=file.path(dir,'assembly_rows_lengths.csv'), sep=',', row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(! is.na(result$disassembly$length))
	if (! is.null(row_nums)) {
		rows_and_length = cbind(row_nums, result$disassembly$length[row_nums]);
		write.table(rows_and_length,file=file.path(dir,'disassembly_rows_lengths.csv'), sep=',', row.names=FALSE, col.names=FALSE)
	}
}

write_high_r_rows <- function(result, dir, file=c('assembly_R_sq.csv','disassembly_R_sq.csv','neg_slope_R_sq.csv'), min_R_sq = 0.9) {
	if (! file.exists(dir)) {
		dir.create(dir,recursive=TRUE)
	}
	
	row_nums = which(is.finite(result$assembly$R_sq) & result$assembly$R_sq > min_R_sq)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[1]), row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(is.finite(result$assembly$R_sq) & result$disassembly$R_sq > min_R_sq)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[2]), row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(is.finite(result$assembly$R_sq) & result$assembly$R_sq > min_R_sq & result$assemblyassemblyassembly$slope < 0)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[3]), row.names=FALSE, col.names=FALSE)
	}
}

write_model_data_to_csv <- function(model, output_file, pre_filt = FALSE) {
    props_to_keep = c('slope','R_sq','p_val','length');
    compiled_props = list();
    for (type in c('assembly','disassembly')) {
        model_set = subset(model[[type]],!is.na(R_sq))
        model_set$lineage_nums = c(model_set$lineage_nums, which(!is.na(model[[type]]$R_sq)));
        if (pre_filt) {
            model_set = subset(model_set, p_val < 0.05 & slope > 0);
        }
        compiled_props$type = c(compiled_props$type, rep(type, dim(model_set)[[1]]));
        compiled_props$lineage_num = c(compiled_props$lineage_num, model_set$lineage_num);
        for (i in props_to_keep) {
            compiled_props[[i]] = c(compiled_props[[i]], model_set[[i]]);
        }
    }
    compiled_props = as.data.frame(compiled_props);
    write.csv(compiled_props, file=output_file, row.names=F)
}
########################################
#Spacial Functions
########################################
correlate_signal_vs_dist <- function(intensities, cent_x, cent_y, min_overlap=40) {
	correlations = c()
	distances = c()
	for (i in 1:dim(intensity_data)[[1]]) {		
		this_row = intensity_data[i,];
		if (sum(! is.na(this_row)) < min_overlap) {
			next
		}
		overlap_rows = which(apply(intensity_data,1, function(data_2) enough_overlap(this_row, data_2, min_overlap=min_overlap)));		
		overlap_rows = setdiff(overlap_rows,1:i)
		
		for (j in overlap_rows) {
			overlap_entries = which(! is.na(this_row) & ! is.na(intensity_data[j,]))
			correlations = c(correlations, cor(as.numeric(this_row[overlap_entries]), as.numeric(intensity_data[j,overlap_entries])))
			
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
            write_model_data_to_csv(average_model[[1]], file.path(data_dir,'..','assemb_disassem_rates.csv'))
            write_model_data_to_csv(average_model[[1]], file.path(data_dir,'..','assemb_disassem_rates_filtered.csv'))
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
