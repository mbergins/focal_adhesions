###############################################################################
# Assembly/Disassembly Rate models
###############################################################################

build_bilinear_models <- function(data,exp_props,min.phase.length=10,time.spacing=1,
    diagnostic.figure = NA) {
    
    models = list()
    
    FA_longevity = rowSums(is.finite(data));
    meets_longev = FA_longevity >= 2*min.phase.length;

    birth_filter = is.nan(data[,1]) & ! as.logical(exp_props$split_birth_status) & meets_longev; 
    
    death_filter = is.nan(data[,dim(data)[2]]) & as.logical(exp_props$death_status) & meets_longev; 
    
    stability_filter = birth_filter & death_filter;
    
    lineages_to_analyze = which(birth_filter | death_filter);
    
    for (FA_number in lineages_to_analyze) {
        #find the models where the birth or death filter properties are met
		
		only.data = na.omit(data[FA_number,]);

        if (birth_filter[FA_number]) {
            time.series = list(value=only.data,
                time=seq(0,along.with=only.data,by=time.spacing));

            assembly_model_results = fit_all_possible_log_models(time.series, 
                min.phase.length=min.phase.length);
        }
        if (death_filter[FA_number]) {
            only.data.reverse = rev(only.data);
            time.series.reverse = list(value=only.data.reverse,
                time=seq(0,along.with=only.data.reverse,by=time.spacing));

            disassembly_model_results = fit_all_possible_log_models(time.series.reverse, 
                min.phase.length=min.phase.length);
        }
        
        #find the best models
        if (birth_filter[FA_number] && death_filter[FA_number]) {
            best_indexes = find_optimum_model_indexes(assembly_model_results,disassembly_model_results);
            
            best_assembly_model = assembly_model_results[best_indexes[1],];
            best_assembly_model$FA_number = FA_number;
            models$assembly = rbind(models$assembly, best_assembly_model);
            
            best_disassembly_model = disassembly_model_results[best_indexes[2],];
            best_disassembly_model$FA_number = FA_number;
            models$disassembly = rbind(models$disassembly, best_disassembly_model);
            
            if (! is.na(best_assembly_model$image_count) & 
                ! is.na(best_disassembly_model$image_count)) {
              stability_model = find_stability_properties(time.series,
                                                          best_assembly_model,
                                                          best_disassembly_model);
              stability_model$FA_number = FA_number;
              stability_model$half_peak_time = find_time_between_half_peaks(only.data,
                                                                            stability_model$mean_intensity,
                                                                            assemb_model=best_assembly_model,
                                                                            disassem_model=best_disassembly_model);
              
              models$stability = rbind(models$stability, stability_model);
            }
        } else if (birth_filter[FA_number]) {
            best_indexes = find_optimum_model_indexes(assembly=assembly_model_results,disassembly=NULL);
            
            best_assembly_model = assembly_model_results[best_indexes[1],];
            best_assembly_model$FA_number = FA_number;
            models$assembly = rbind(models$assembly, best_assembly_model);
        } else if (death_filter[FA_number]) {
            best_indexes = find_optimum_model_indexes(assembly=NULL,disassembly=disassembly_model_results);
            
            best_disassembly_model = disassembly_model_results[best_indexes[2],];
            best_disassembly_model$FA_number = FA_number;
            models$disassembly = rbind(models$disassembly, best_disassembly_model);
        }
    }
    
    models$assembly = as.data.frame(models$assembly);
    models$disassembly = as.data.frame(models$disassembly);
    models$stability = as.data.frame(models$stability);
    models$exp_data = data;
    models$exp_props = exp_props;
    
    print(paste('Built ', length(which(!is.na(models$assembly$slope))), 
        ' assembly models, built ', length(which(!is.na(models$disassembly$slope))), 
        ' disassembly models',sep=''))

    return(models)
}

fit_all_possible_log_models <- function(time.series,min.phase.length=10) {
    data_summary = list();
    
    data_summary$min.val = min(time.series$value,na.rm=T);
    data_summary$max.val = max(time.series$value,na.rm=T);

    time.series$value = time.series$value/time.series$value[1];
    time.series$value[time.series$value < 0] = 1E-5;

    for (i in min.phase.length:length(time.series$value)) {
        this_model = lm(log(time.series$value[1:i]) ~ time.series$time[1:i]);
        model_summary = summary(this_model);
        
        data_summary$adj.r.squared = c(data_summary$adj.r.squared, model_summary$adj.r.squared);
        data_summary$r.squared = c(data_summary$r.squared, model_summary$r.squared);
        data_summary$p.value = c(data_summary$p.value, model_summary$coefficients[2,4]);
        data_summary$slope = c(data_summary$slope, coef(this_model)[2]);
        data_summary$time_length = c(data_summary$time_length, time.series$time[i]);
        data_summary$image_count = c(data_summary$image_count, i);
    }

    data_summary = as.data.frame(data_summary);
}

find_time_between_half_peaks <- function(fa_intensity,stability_mean,
										 assemb_model,disassem_model) {

	#Requested by Pradeep Uchill to implement a metric from a Maddy Parson
	#paper.
	fa_intensity_min = fa_intensity - min(fa_intensity);
	
	stability_mean_min = stability_mean - min(fa_intensity);
	
	meet_half_stability = which(fa_intensity_min > stability_mean_min*0.5)

	time_between = max(meet_half_stability) - min(meet_half_stability);
	return(time_between);
}

find_optimum_model_indexes <- function(assembly=NULL,disassembly=NULL) {
    best_fit_indexes = c(NA, NA);
   
    if (is.null(disassembly)) {
        best_index = min(which(assembly$adj.r.squared == max(assembly$adj.r.squared)))
        best_fit_indexes[1] = best_index
    } else if (is.null(assembly)) {
        best_index = min(which(disassembly$adj.r.squared == max(disassembly$adj.r.squared)));
        best_fit_indexes[2] = best_index
    } else {
        r_sums = matrix(NA,length(assembly$slope), length(disassembly$slope));
        
        best_indexes = c();
        
        max_image_count = max(assembly$image_count);

        for (as_index in 1:dim(r_sums)[1]) {
            for (dis_index in 1:dim(r_sums)[2]) {
                #skip calculating R sum if the total length of the combination
                #is greater than the length of the time series set
                total_images_in_models = assembly$image_count[as_index] + 
                                         disassembly$image_count[dis_index];
                if (total_images_in_models > max_image_count) {
                    next;
                }

                r_sums[as_index,dis_index] = assembly$adj.r.squared[as_index] + disassembly$adj.r.squared[dis_index]
            }
        }
        
        best_fit_indexes = which(max(r_sums,na.rm=T) == r_sums, arr.ind=T);
        
        #dealing with the unlikely case that there is more than one hit on the
        #highest summed R-squared value, fortunately, "which" returns the set
        #of indexes closest to the upper right corner first, that happens to be
        #the shortest model set and also the one we want
        if (dim(best_fit_indexes)[1] > 1) {
            best_fit_indexes = best_fit_indexes[1,];
        }
    }

    return(best_fit_indexes);
}

find_stability_properties <- function(time.series,best_assembly,best_disassembly) {
    stability_props = list();
    stability_props$image_count = length(time.series$value) - 
        best_assembly$image_count - 
        best_disassembly$image_count;
    
	#If the stability phase contains at least one image, we will only use those
	#images to make the remaining calculation, otherwise, use the last data
	#point in assembly and the first in disassembly
    if (stability_props$image_count >= 1) {
        stability_indexes = seq(best_assembly$image_count+1, length=stability_props$image_count);
        stability_values = time.series$value[stability_indexes];
	} else {
        stability_indexes = c(best_assembly$image_count, best_assembly$image_count+1);
        stability_values = time.series$value[stability_indexes];
	}

	stability_fold_change = stability_values/mean(c(time.series$value[1],tail(time.series$value,n=1)));
	stability_props$mean_intensity = mean(stability_values);
	stability_props$mean_fold_change = mean(stability_fold_change);
	stability_props$stdev = sd(stability_fold_change);
	stability_props$coeff_of_var = stability_props$stdev/stability_props$mean_fold_change;

    stability_props = as.data.frame(stability_props);
    return(stability_props);
}

###########################################################
# Plotting
###########################################################

draw_diagnostic_traces <- function(models,file.name,time.spacing=1) {
    ad_nums = sort(union(models$assembly$FA_number,models$disassembly$FA_number));
    
    y_limits = c(min(models$assembly$min,models$disassembly$min,na.rm=T),
        max(models$assembly$max,models$disassembly$max,na.rm=T));

    pdf(file.name)
    for (this_FA_num in ad_nums) {
        phase_lengths = c(NA,NA);
        R_sq = c(NA,NA);
        slopes = c(NA,NA);

        if (any(this_FA_num == models$assembly$FA_number)) {
            this_assemb = subset(models$assembly, FA_number == this_FA_num);
            phase_lengths[1] = this_assemb$image_count;
            R_sq[1] = this_assemb$r.squared;
            slopes[1] = this_assemb$slope;
        }
        if (any(this_FA_num == models$disassembly$FA_number)) {
            this_disassemb = subset(models$disassembly, FA_number == this_FA_num);
            phase_lengths[2] = this_disassemb$image_count;
            R_sq[2] = this_disassemb$r.squared;
            slopes[2] = this_disassemb$slope;
        }
            
        plot_ad_intensity(models,this_FA_num,phase_lengths,R_sq,y_limits,slopes,time.spacing);
    }
    graphics.off()
}

produce_rate_filters <- function(single.exp.data, model_count, min.r.sq=-Inf, 
    max.p.val=0.05, pos.slope=TRUE, old.names=F, min.cent.dist=NA) {
    
    filter_sets = list()
    
    # First we cycle through all the variables used to filter the assembly and
    # disassembly phases, collecting logic table describing which adhesions
    # fulfill which criteria
    vars_to_filter = c("assembly", "disassembly")
    for (phase in vars_to_filter) {
        #apply a bonferroni correction to the maximum acceptable p-value
        #calculation if model count provided
        corrected.p.value = max.p.val;
        if (! is.na(model_count)) {
            corrected.p.value = corrected.p.value / model_count;
        }
        
		filter_sets[[phase]]$good.r.sq = ! is.na(single.exp.data[[phase]]$adj.r.squared) & 
			single.exp.data[[phase]]$adj.r.squared >= min.r.sq;
		filter_sets[[phase]]$low.p.val = ! is.na(single.exp.data[[phase]]$p.value) & 
			single.exp.data[[phase]]$p.value < corrected.p.value; 
		filter_sets[[phase]]$pos.slope = ! is.na(single.exp.data[[phase]]$slope) & 
			single.exp.data[[phase]]$slope > 0;

        if (! is.na(min.cent.dist)) {
            thresh = quantile(single.exp.data$exp_props$Mean_FA_cent_dist,c(min.cent.dist))
            filter_sets[[phase]]$cent_dist = ! is.na(single.exp.data$exp_props$Mean_FA_cent_dist) & 
                single.exp.data$exp_props$Mean_FA_cent_dist > thresh; 
        }
    }
    
    # Now all the filters are cascaded to produce the final filter set, then
    # the pos.slope variable is checked and applied if requested
    final_filters = list()
    for (phase in vars_to_filter) {
        final_filters[[phase]] = (filter_sets[[phase]]$good.r.sq & filter_sets[[phase]]$low.p.val);
        if (any(names(filter_sets[[phase]]) == "cent_dist")) {
            final_filters[[phase]] = (final_filters[[phase]] & filter_sets[[phase]]$cent_dist);
        }
    }
    
	assemb_FA_nums = single.exp.data$assem$FA_number[final_filters$assembly];
	dis_FA_nums = single.exp.data$dis$FA_number[final_filters$disassembly];
    joint_FA_nums = intersect(assemb_FA_nums,dis_FA_nums);
	final_filters$joint = is.element(single.exp.data$stability$FA_num,joint_FA_nums)
    
	for (filter_type in names(final_filters)) {
        final_filters[[filter_type]] = as.logical(final_filters[[filter_type]])
    }

    return(final_filters)
}

find_number_of_models <- function(model_set) {
    count = 0;
    for (i in 1:length(model_set)) {
        count = count + sum(! is.na(model_set[[i]]$assembly$p.value));
        count = count + sum(! is.na(model_set[[i]]$disassembly$p.value));
    }
    return(count)
}

write_assembly_disassembly_periods <- function(ad_kinetics, dir) {
	if (! file.exists(dir)) {
		dir.create(dir,recursive=TRUE)
	}
	
    just_assemb = subset(ad_kinetics,class=="Assembly");
    rows_and_length = cbind(just_assemb$FA_number,just_assemb$image_count);
    write.table(rows_and_length,file=file.path(dir,'assembly_rows_lengths.csv'), 
                sep=',', row.names=FALSE, col.names=FALSE)
    
    just_dis = subset(ad_kinetics,class=="Disassembly");
    rows_and_length = cbind(just_dis$FA_number,just_dis$image_count);
    write.table(rows_and_length,file=file.path(dir,'disassembly_rows_lengths.csv'), 
                sep=',', row.names=FALSE, col.names=FALSE)
}

produce_simple_CSV_output_set <- function(models,time_spacing) {
    output_properties = c('class','FA_number','slope','p.value','r.squared',
                          'adj.r.squared','image_count');

    #Assembly Models Extraction
    assembly_models = models$assembly;
	if (dim(assembly_models)[1] > 0) {
		assembly_models$class = rep('Assembly',dim(assembly_models)[1]);
		assembly_adhesions = subset(assembly_models, 
									!is.na(p.value) & p.value < 0.05 & slope > 0, 
									select = output_properties);
		assembly_adhesions$slope = round(assembly_adhesions$slope,4);
		assembly_adhesions$adj.r.squared = round(assembly_adhesions$adj.r.squared,3);
		# assembly_adhesions$phase_length = assembly_adhesions$image_count * time_spacing;
		# assembly_adhesions$image_count <- NULL;
	} else {
		assembly_adhesions = NULL;
	}

    #Disassembly Models Extraction
	disassembly_models = models$disassembly;
	if (dim(disassembly_models)[1] > 0) {
		disassembly_models$class = rep('Disassembly',dim(disassembly_models)[1]);
		disassembly_adhesions = subset(disassembly_models, 
									   !is.na(p.value) & p.value < 0.05 & slope > 0, 
									   select = output_properties);
		disassembly_adhesions$slope = round(disassembly_adhesions$slope,4);
		disassembly_adhesions$adj.r.squared = round(disassembly_adhesions$adj.r.squared,3);
		# disassembly_adhesions$phase_length = disassembly_adhesions$image_count * time_spacing;
		# disassembly_adhesions$image_count <- NULL;
	} else {
		disassembly_models = NULL;
	}

    #Stability Data Set Extraction
    valid_stability_nums = intersect(assembly_adhesions$FA_number,disassembly_adhesions$FA_number);
    stability_adhesions = models$stability[match(valid_stability_nums,models$stability$FA_number),];
    stability_adhesions_full_data = stability_adhesions;

    stability_adhesions$mean_intensity <- NULL;
    stability_adhesions$mean_fold_change <- NULL;
    stability_adhesions$class = rep('Stability',length(valid_stability_nums));

    for (prop in output_properties) {
        if (! any(names(stability_adhesions) == prop)) {
            stability_adhesions[[prop]] = rep(NA,length(valid_stability_nums));
        }
    }
    
    stability_adhesions = stability_adhesions[names(assembly_adhesions)];

    output_sets = list();
    output_sets$ad_kinetics = rbind(assembly_adhesions,disassembly_adhesions,stability_adhesions);
    output_sets$ad_kinetics$phase_length = output_sets$ad_kinetics$image_count * time_spacing;
	
	stability_adhesions_full_data$half_peak_time = stability_adhesions_full_data$half_peak_time * time_spacing;
	output_sets$stability_props = stability_adhesions_full_data;
    
    return(output_sets);
}

################################################################################
# Main Command Line Program
################################################################################

args = commandArgs(TRUE);
if (length(args) != 0) {
    debug = FALSE;
    min_length = 10;
    time_spacing = 1;
    
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
    class(time_spacing) <- "numeric";
    if (exists('data_dir') & exists('data_file')) {
        #######################################################################
        # Reading in data files
        #######################################################################
        if (file.exists(file.path(data_dir,'lin_time_series',data_file))) {
            data_set = as.matrix(read.csv(file.path(data_dir,'lin_time_series',data_file),header=F));
        } else {
            print(paste('Could not find data file: ', file.path(data_dir,'lin_time_series',data_file)))
            stop()
        }

        if (file.exists(file.path(data_dir,'single_lin.csv'))) {
            exp_props = read.csv(file.path(data_dir,'single_lin.csv'));
        } else {
            print(paste('Could not find exp props file: ', file.path(data_dir,'single_lin.csv')))
            stop()
        }
        
        if (file.exists(file.path(data_dir,'lin_time_series','FA_angle_recentered.csv'))) {
            angle_recenter = as.matrix(read.csv(file.path(data_dir,'lin_time_series','FA_angle_recentered.csv'),header=F));
            exp_props$Mean_FA_recentered_angle = rowMeans(angle_recenter,na.rm=T);
        } else {
            print(paste('Could not find FA_angle_recentered.csv'))
        }
        
        output_folder = file.path(data_dir,'assem_disassem_models');
		dir.create(output_folder,recursive=TRUE,showWarnings=F);
        
        #######################################################################
        # Model Building and Output
        #######################################################################
        models = build_bilinear_models(data_set,exp_props, min.phase.length = min_length, 
            time.spacing = time_spacing);
        models$exp_dir = data_dir;

		R_models_file = sub("\\..*",".Rdata",data_file,perl=T)
        output_file = file.path(output_folder,R_models_file);
        save(models,file = output_file)
        
        diagnostic_diagrams_file = sub("\\..*",".pdf",data_file,perl=T)
		output_file = file.path(output_folder,diagnostic_diagrams_file);
        source('FA_analysis_lib.R')
        draw_diagnostic_traces(models,output_file,time.spacing=time_spacing);
        
        #######################################################################
        # Simple CSV File Output
        #######################################################################
        output_sets = produce_simple_CSV_output_set(models,time_spacing)
        for (output_type in names(output_sets)) {
            write.csv(output_sets[[output_type]],
                file=file.path(data_dir,paste0(output_type,'.csv')),row.names=F)
        }

        write_assembly_disassembly_periods(output_sets$ad_kinetics,output_folder);
    }
}
