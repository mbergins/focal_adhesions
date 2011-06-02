###############################################################################
# FA Orientation Models
###############################################################################

gather_FA_orientation_data <- function(exp_dir,fixed_best_angle = NA,
    min.ratio = 3,output_file = 'FA_orientation.Rdata',diagnostic.figure=F) {

    data_set = read_in_orientation_data(exp_dir, min.ratio=min.ratio);
    data_set$angle_search = test_dom_angles(data_set$subseted_data$orientation);

    if (is.na(fixed_best_angle)) {
        data_set$best_angle = find_best_alignment_angle(data_set$angle_search)
    } else {
        data_set$best_angle = fixed_best_angle
        data_set$actual_best_angle = find_best_alignment_angle(data_set$angle_search)
    }
    data_set$corrected_orientation = apply_new_orientation(data_set$subseted_data$orientation,
        data_set$best_angle)
    
    per_image_dom_angle = find_per_image_dom_angle(data_set$mat, min.ratio=min.ratio)
    data_set$per_image_dom_angle = per_image_dom_angle

    data_set$single_ads = analyze_single_adhesions(data_set)
    
    save(data_set,file=file.path(exp_dir,'..',output_file))
    
    if (! diagnostic.figure) {
        return(data_set);
    }

    #diagnostic figure
    pdf(file.path(exp_dir,'..','adhesion_orientation.pdf'), height=7*(3/2))
    layout(rbind(c(1,2),c(3,4),c(5,6)))
    par(bty='n', mar=c(4,4.2,2,0),mgp=c(2,1,0))
    
    hist(data_set$subseted_data$orientation,main='Pos X-axis Reference',
        xlab=paste('Angle n=',dim(data_set$subseted_data)[1],sep=''), 
		xlim=c(-100,100));
    
    plot(data_set$angle_search$test_angles,data_set$angle_search$angle_FAAI,typ='l',
        xlab='Dominant Search Angle',ylab='FA Alignment Index',ylim=c(0,90));
    lines(data_set$angle_search$test_angles, abs(data_set$angle_search$mean_angle), col='red')
    lines(data_set$angle_search$test_angles, abs(data_set$angle_search$median_angle), col='blue')
    if (! is.na(fixed_best_angle)) {
        segments(fixed_best_angle,0,fixed_best_angle,2000,col='blue')
        segments(data_set$actual_best_angle,0,data_set$actual_best_angle,2000,col='green')
    } else {
        segments(data_set$best_angle,0,data_set$best_angle,2000,col='green')
    }

    hist(data_set$corrected_orientation,
        main=paste('Rotated ',data_set$best_angle,'\u00B0 / ',
            sprintf('%0.1f',90-sd(data_set$corrected_orientation)),' FAAI',sep=''),
        xlab=paste('Angle n=',dim(data_set$subseted_data)[1],sep=''),
		xlim=c(-100,100));
    
    plot(per_image_dom_angle,xlab='Image Number',ylab='Dominant Angle',ylim=c(0,180));

    hist(data_set$single_ads$best_angle,main='Single Adhesions - Best Angle',
        xlab=paste('Angle n=',length(data_set$single_ads$best_angle),sep=''),xlim=c(0,180));
    
    hist(data_set$single_ads$FAAI,main='Single Adhesions - Max FAAI',
        xlab=paste('Angle n=',length(data_set$single_ads$best_angle),sep=''),xlim=c(0,90));

    graphics.off()

    return(data_set)
}

read_in_orientation_data <- function(time_series_dir,min.ratio = 3) {
    data_set = list();
    data_set$mat$orientation = read.csv(file.path(time_series_dir,'Orientation.csv'),header=F);
    data_set$mat$area = read.csv(file.path(time_series_dir,'Area.csv'),header=F);
    major_axis = read.csv(file.path(time_series_dir,'MajorAxisLength.csv'),header=F);
    minor_axis = read.csv(file.path(time_series_dir,'MinorAxisLength.csv'),header=F);
    
    data_set$mat$ratio = major_axis/minor_axis;

    unlist_data_set = list()
    for (i in names(data_set$mat)) {
        unlist_data_set[[i]] = unlist(data_set$mat[[i]]);
    }
    unlist_data_set = as.data.frame(unlist_data_set);
    data_set$subseted_data = subset(unlist_data_set, !is.nan(unlist_data_set$ratio) & 
        unlist_data_set$ratio >= min.ratio);
    
    return(data_set);
}

find_per_image_dom_angle <- function(mat_data, min.ratio=3) {
    best_angles = c()
    for (i_num in 1:dim(mat_data$orientation)[2]) {
        this_orientation = mat_data$orientation[,i_num];
        this_ratio = mat_data$ratio[,i_num];

        good_rows = !is.na(this_orientation) & this_ratio >= min.ratio;
        
        angle_search = test_dom_angles(this_orientation[good_rows]);
        best_angle = find_best_alignment_angle(angle_search)

        best_angles = c(best_angles, best_angle);
    }

    return(best_angles)
}

test_dom_angles <- function(orientation, search_resolution = 0.1) {
    # I'll include and then discard the last angle because 180 is the same as 0
    # degrees rotation in this case, but we need a consistant end point to allow
    # the search resolution to vary
    angles_to_test = seq(0,180,by=search_resolution)
    angles_to_test = angles_to_test[-length(angles_to_test)];
    
    results = list(x = angles_to_test, test_angles = angles_to_test);

    new_angle_FAAI = c()
    mean_angle = c()
    median_angle = c()
    for (angle in angles_to_test) {
        new_orientation = apply_new_orientation(orientation,angle);
        mean_angle = c(mean_angle,mean(new_orientation, na.rm=T));
        median_angle = c(median_angle,median(new_orientation, na.rm=T));
        new_angle_FAAI = c(new_angle_FAAI, find_FAAI_from_orientation(new_orientation));
    }
        
    results$y = new_angle_FAAI;
    results$mean_angle = mean_angle;
    results$median_angle = median_angle;
    results$angle_FAAI = new_angle_FAAI;

    return(results)
}

find_FAAI_from_orientation <- function(orientation_data) {
    return(90-sd(orientation_data, na.rm=T));
}

find_best_alignment_angle <- function(test_angle_set) {
    # to select the best angle, find the angles with the maximum FAAI, from
    # those angles select the angle with the lowest (meaning closest to zero)
    # mean angle
    max_FAAI = max(test_angle_set$angle_FAAI)
    max_FAAI_indexes = which(max_FAAI == test_angle_set$angle_FAAI)
    
    abs_mean = abs(test_angle_set$mean_angle)[max_FAAI_indexes];
    sorted_abs_mean = sort(abs_mean,decreasing = F, index.return = T)
    
    best_index = max_FAAI_indexes[sorted_abs_mean$ix[1]]

    best_angle = test_angle_set$test_angles[best_index];
    return(best_angle)
}

apply_new_orientation <- function(orientation_data,angle) {
    orientation_data = orientation_data - angle;
    
    less_neg_ninety = ! is.na(orientation_data) & orientation_data < -90;
    orientation_data[less_neg_ninety] = orientation_data[less_neg_ninety] + 180;

    return(orientation_data)
}

analyze_single_adhesions <- function(align_data, min.data.points = 5, min.ratio = 3) {
    orientations = as.matrix(align_data$mat$orientation);
    ratio = as.matrix(align_data$mat$ratio);

    num_above_ratio_limit = apply(ratio,1,function(x) sum(! is.na(x) & x >= min.ratio));

    passed_ad_nums = which(num_above_ratio_limit >= min.data.points);
    
    single_ad_data = list()
    single_ad_data$ad_nums = passed_ad_nums

    time = 0:(dim(orientations)[2]-1)
    for (ad_num in passed_ad_nums) {
        this_ad_filter = ! is.na(ratio[ad_num,]) & ratio[ad_num,] >= min.ratio;

        these_angles = orientations[ad_num,this_ad_filter];
        these_ratio = ratio[ad_num,this_ad_filter];
        stopifnot(length(these_angles) == num_above_ratio_limit[ad_num])
        
        angle_search = test_dom_angles(these_angles);
        best_angle = find_best_alignment_angle(angle_search);

        single_ad_data$best_angle = c(single_ad_data$best_angle,best_angle)
        single_ad_data$FAAI = c(single_ad_data$FAAI, max(angle_search$angle_FAAI))
        
        lin_model = lm(orientations[ad_num,]~time);
        lin_model_summary = summary(lin_model);

        single_ad_data$lin_model_p = c(single_ad_data$lin_model_p,lin_model_summary$coefficients[2,4]);
    }

    return(single_ad_data);
}

###########################################################
# Processing Single Adhesion Data
###########################################################

filter_alignment_data <- function(align_data, min.data.points = 10, min.ratio = 3, 
	min.area = 60) {
    
    orientation = as.matrix(align_data$mat$orientation);
    ratio = as.matrix(align_data$mat$ratio);
    area = as.matrix(align_data$mat$area);
    
    above_ratio_limit = ! is.na(ratio) & ratio >= min.ratio
    above_area_limit = ! is.na(area) & area >= min.area

    above_all_limits = above_ratio_limit & above_area_limit;
    num_above_limits = rowSums(above_all_limits)
	# conse_above_limits = apply(above_all_limits,1,number_consecutive_trues)
    passed_ad_nums = which(num_above_limits >= min.data.points)

    if (any(names(align_data) == "lineage_data")) {
		if (any(names(align_data$lineage_data) == "split_count")) {
			passed_ad_nums = intersect(passed_ad_nums, which(align_data$lineage_data$split_count <= 5));
		}
		if (any(names(align_data$lineage_data) == "merge_count")) {
			passed_ad_nums = intersect(passed_ad_nums, which(align_data$lineage_data$merge_count <= 5));
		}
    }

	for (ad_num in 1:dim(orientation)[1]) {
		if (any(passed_ad_nums == ad_num)) {
			next;
		} else {
			orientation[ad_num,] = NA;
		}
	}

	orientation[! above_all_limits] = NA;
	
    align_data$mat$filtered_orientation = orientation;

    return(align_data);
}

plot_single_adhesion_orientations <- function(align_data, min.data.points = 20,
 	max.data.points=Inf, min.ratio = 3, dominant.angle = 0, which.ads=NA, 
 	min.area = 100,main=NA,with.tags=F) {

    orientations = as.matrix(align_data$mat$orientation);
    ratio = as.matrix(align_data$mat$ratio);
    area = as.matrix(align_data$mat$area);
 
    align_data = build_single_orientation_filters(align_data,
        min.data.points = min.data.points, max.data.points=max.data.points,
        min.ratio = min.ratio, min.area = min.area)
	
	passed_ad_nums = which(rowSums(! is.na(aligned_data$mat$filtered_orientations)) != 0)
	
    if (! is.na(which.ads)) {
        passed_ad_nums = which.ads
    }
    
    if (length(passed_ad_nums) == 0) {
        return();
    }

    colors = rainbow(length(passed_ad_nums));
    colors = sample(colors,length(colors),replace=F);

    time = (0:(dim(orientations)[2]-1))*2.5
    count = 1;
    for (ad_num in passed_ad_nums) {
        this_ad_filter = above_all_limits[ad_num,];
        
        angle_sequence = orientations[ad_num,];
        angle_sequence[! this_ad_filter] = NA;
        angle_sequence = apply_new_orientation(angle_sequence,dominant.angle);
        
        if (count == 1) {
            xlims = c(0,max(time))
            if (with.tags) {
                xlims = c(0,max(time)*1.1);
            }

            plot(time,angle_sequence,ylim=c(-90,90),col=colors[count],
                pch=19,cex=0.25,xlim=xlims,axes=F,
                xlab='Time (min)',ylab='FA Orientation (degrees)');
            lines(time,angle_sequence, col=colors[count],pch=19,cex=0.25)

            if (! is.na(main)) {
                title(main=main);
            }
            
            p_limits = par("usr");
            text(p_limits[2],p_limits[4],
                paste(align_data$best_angle,'\u00B0','/','n=',length(passed_ad_nums),sep=''),
                pos=2,offset=c(0,-1));

            axis(1)
            axis(2, at = seq(-90,90,by=90))
        } else {
            points(time,angle_sequence, col=colors[count],pch=19,cex=0.25)
            lines(time,angle_sequence, col=colors[count],pch=19,cex=0.25)
        }

        ad_present_time = (which(!is.na(angle_sequence))-1)*2.5;
        omited_sequence = na.omit(angle_sequence);
        # segments(min(ad_present_time),mean(angle_sequence,na.rm=T),
        #     max(ad_present_time),mean(angle_sequence,na.rm=T),
        #     col=colors[count])

        lin_model = lm(angle_sequence~time);
        lin_model_summary = summary(lin_model);
        # lines(ad_present_time,predict(lin_model))
        
        if (with.tags) {
            single_ad_index = which(align_data$single_ads$ad_nums == ad_num);
            # text(max(ad_present_time),mean(angle_sequence,na.rm=T),
            #     sprintf('%.1f',align_data$single_ads$FAAI[single_ad_index]),col=colors[count],
            #     pos=3)
            text(max(ad_present_time),tail(omited_sequence,1),
                ad_num,col=colors[count],
                pos=1)
            # text(max(ad_present_time),mean(angle_sequence,na.rm=T),
            #     sprintf('%.1e',lin_model_summary$coefficients[2,4]),col=colors[count],
            #     pos=4)
        }
        count = count + 1;
    }
    pl_size = par("usr");
    segments(pl_size[1],0,pl_size[2],0)

    print(count)
}

adhesion_angle_deviance <- function(orientations,min.data.points) {
	passed_ads = which(rowSums(! is.na(orientations)) >= 1)
	mean_dev = c()
	for (ad_num in passed_ads) {
		or_set = orientations[ad_num,];
		or_set = na.omit(or_set);
        
		angle_test = test_dom_angles(or_set);
		best_angle = find_best_alignment_angle(angle_test);

		or_set = apply_new_orientation(or_set, best_angle);

		diffs = abs(or_set[2:length(or_set)] - or_set[1]);
		mean_dev = c(mean_dev,mean(diffs));
	}
	return(data.frame(mean_dev = mean_dev, ad_num=passed_ads))
}

number_consecutive_trues <- function(logical_seq) {
    logical_seq = as.logical(logical_seq)

    max_consec = 0
    this_consec = 0
    for (i in 1:length(logical_seq)) {
        if (logical_seq[i]) {
            this_consec = this_consec + 1;
        } else {
            if (this_consec > max_consec) {
                max_consec = this_consec;
            } else {
                this_consec = 0;
            }
        }
    }
    
    #if the entire sequence is true, the else clause above won't be hit, check
    #for this and set max_consec appropriately
    if (this_consec > max_consec) {
        max_consec = this_consec
    }

    return(max_consec);
}

stopifnot(number_consecutive_trues(c(rep(F,10),rep(T,10),rep(F,50),rep(T,200))) == 200)
stopifnot(number_consecutive_trues(c(rep(F,10),rep(T,10),rep(F,50),rep(T,20))) == 20)
stopifnot(number_consecutive_trues(c(rep(T,25),rep(F,10),rep(T,10),rep(F,50),rep(T,20))) == 25)
stopifnot(number_consecutive_trues(rep(T,20)) == 20)

gather_all_single_adhesion_deviances <- function(align_files, lin_files, min.area=-Inf, min.data.points=2) {
	all_dev_data = list();
	for (i in 1:length(align_files)) {
		sample_data = get(load(align_files[i]))
		sample_data$lineage_data <- read.table(lin_files[i],sep=',',header=T);
		
		sample_data_filtered = filter_alignment_data(sample_data, 
			min.area=min.area, min.data.points=min.data.points);
		overall_dev = adhesion_angle_deviance(sample_data_filtered$mat$filtered_orientation);
		overall_dev$exp_num = rep(i,dim(overall_dev)[1])
		overall_dev$longevity = sample_data$lineage_data$longevity[overall_dev$ad_num];

		all_dev_data = rbind(all_dev_data, overall_dev);
        print(paste('Done with ',i,'/',length(align_files),sep=''));
	}
    return(all_dev_data)
}

get_mean_major_minor_ratio <- function(lin_ts_folder,min.longevity = 10) {
    major_file = Sys.glob(file.path(lin_ts_folder, 'Major*'))
    minor_file = Sys.glob(file.path(lin_ts_folder, 'Minor*'))

    exp_props = read.csv(file.path(lin_ts_folder, '../single_lin.csv'))
    filt = ! exp_props$split_birth_status & exp_props$death_status

    major = as.matrix(read.csv(major_file,header=F));
    minor = as.matrix(read.csv(minor_file,header=F));

    ratio = major / minor

    longevity = rowSums(! is.na(ratio));

    ratio = ratio[longevity > min.longevity & filt,];
    
    return(rowMeans(ratio,na.rm=T))
}

get_mean_eccen <- function(lin_ts_folder,min.longevity = 10) {
    major_file = Sys.glob(file.path(lin_ts_folder, 'Major*'))
    minor_file = Sys.glob(file.path(lin_ts_folder, 'Minor*'))

    exp_props = read.csv(file.path(lin_ts_folder, '../single_lin.csv'))
    filt = ! exp_props$split_birth_status & exp_props$death_status

    major = as.matrix(read.csv(major_file,header=F));
    minor = as.matrix(read.csv(minor_file,header=F));
    
    ratio = (minor/2)/(major/2)
    
    eccen = sqrt(1-ratio^2)

    longevity = rowSums(! is.na(eccen));

    eccen = eccen[longevity > min.longevity & filt,];
    # eccen = eccen[longevity > min.longevity,];
    
    return(rowMeans(eccen,na.rm=T))
}

###########################################################
# Spatial
###########################################################

determine_mean_dist_between <- function(centroid_x,centroid_y,min.overlap=1) {
    dists = matrix(NA,nrow=dim(centroid_x)[1],ncol=dim(centroid_x)[1]);
    overlap_counts = matrix(NA,nrow=dim(centroid_x)[1],ncol=dim(centroid_x)[1]);
    
    ad_present = ! is.na(centroid_x);
    
    for (ad_num in 1:(dim(centroid_x)[1]-1)) {
        for (other_ad_num in (ad_num+1):dim(centroid_x)[1]) {
            overlap_count = sum(ad_present[ad_num,] & ad_present[other_ad_num,]);
            overlap_counts[ad_num,other_ad_num] = overlap_count;
            if (overlap_count >= min.overlap) {
                data_1 = rbind(centroid_x[ad_num,],centroid_y[ad_num,]);
                data_2 = rbind(centroid_x[other_ad_num,],centroid_y[other_ad_num,]);
                
                mean_dist = find_mean_dist(data_1,data_2);
                dists[ad_num,other_ad_num] = mean_dist;
            }
        }
    }
    return(list(dists = dists,overlap_counts = overlap_counts));
}

find_mean_dist <- function(data_1,data_2) {
    mean_dist = NA;
   
    overlap_indexes = !is.na(data_1[1,]) & !is.na(data_2[1,]);
    if (all(! overlap_indexes)) {
        return(mean_dist);
    }

    data_1 = data_1[,overlap_indexes];
    data_2 = data_2[,overlap_indexes];
    
    dists = c()
    for (i in 1:dim(data_1)[2]) {
        dists = c(dists, sqrt((data_1[1,i]-data_2[1,i])^2+(data_1[2,i]-data_2[2,i])^2))
    }
    
    mean_dist = mean(dists);
    return(mean_dist);
}

###########################################################
# Rate and Angle Data
###########################################################

gather_rate_versus_angle_data_set <- function(kinetics_data) {
    align_file = file.path(kinetics_data$exp_dir,'FA_orientation.Rdata');
    
    var_name = load(align_file)
    align_data = get(var_name);
    
    #setting up variables
    single_ads_pos_x_ref = align_data$single_ads$best_angle;
    greater_90 = single_ads_pos_x_ref > 90;
    single_ads_pos_x_ref[greater_90] = 180 - single_ads_pos_x_ref[greater_90]
    
    single_ads_best_ref = apply_new_orientation(single_ads_pos_x_ref,align_data$best_angle);

    temp_birth = kinetics_data$exp_props$birth_i_num;
    temp_birth[is.na(temp_birth)] = 0;
    
    temp_death = kinetics_data$exp_props$death_i_num+1;
    temp_death[is.na(temp_death)] = max(temp_death,na.rm=T);
    longev_unsure = temp_death - temp_birth
    
    #putting the data set together
    combined_data = list()
    ad_nums = align_data$single_ads$ad_nums

    combined_data$ad_nums = ad_nums;
    combined_data$assembly_rate = kinetics_data$assembly$slope[ad_nums]
    combined_data$disassembly_rate = kinetics_data$disassembly$slope[ad_nums]

    combined_data$longevity_unsure = longev_unsure[ad_nums]
    combined_data$longevity = kinetics_data$exp_props$longevity[ad_nums]

    combined_data$angle_dev = single_ads_best_ref
    combined_data$overall_FAAI = rep(max(align_data$angle_search$angle_FAAI),length(ad_nums))
    combined_data$ad_FAAI = align_data$single_ads$FAAI
    
    combined_data = as.data.frame(combined_data)
    
    return(combined_data)
}

get_angle_data_set <- function(kinetics_set) {
    rate_angle_data = list()
    for (i in 1:length(kinetics_set)) {
        temp = gather_rate_versus_angle_data_set(kinetics_set[[i]])
        rate_angle_data = rbind(rate_angle_data, temp)
        # print(paste('Done with ',i,'/',length(kinetics_set),sep=''))
    }
    return(rate_angle_data)
}

split_angle_data <- function(rate_angle_data, angle=45) {
    split_data = list()

    split_data$middle_assembly = subset(rate_angle_data,
        !is.na(assembly_rate) & assembly_rate > 0 & abs(angle_dev) < angle)
    split_data$out_assembly = subset(rate_angle_data,
        !is.na(assembly_rate) & assembly_rate > 0 & abs(angle_dev) >= angle)

    split_data$middle_disassembly = subset(rate_angle_data,
        !is.na(disassembly_rate) & disassembly_rate > 0 & abs(angle_dev) < angle)
    split_data$out_disassembly = subset(rate_angle_data,
        !is.na(disassembly_rate) & disassembly_rate > 0 & abs(angle_dev) >= angle)

    split_data$middle_longev = subset(rate_angle_data, abs(angle_dev) < angle)
    split_data$out_longev = subset(rate_angle_data, abs(angle_dev) >= angle)
    
    return(split_data)
}

###########################################################
# Data Loading
###########################################################

load_alignment_props <- function(alignment_models) {
    if (length(alignment_models) == 0) {
        print('Problem, no alignment models submitted.');
    }

    align_props = list()
    for (align_file in alignment_models) {
        data_set = get(load(align_file));
        
        align_props$best_FAAI = c(align_props$best_FAAI, 
            find_FAAI_from_orientation(data_set$corrected_orientation));
        if (any(names(data_set) == "actual_best_angle")) {
            align_props$actual_best_angle = c(align_props$actual_best_angle, data_set$actual_best_angle);
            align_props$best_angle = c(align_props$best_angle, data_set$best_angle);
        } else {
            align_props$best_angle = c(align_props$best_angle, data_set$best_angle);
        }
        align_props$align_file = c(align_props$align_file, align_file);
        print(paste('Done loading:', align_file))
    }
    align_props = as.data.frame(align_props);
    return(align_props)
}

################################################################################
# Main Command Line Program
################################################################################

args = commandArgs(TRUE);
if (length(args) != 0) {
    debug = FALSE;
    fixed_best_angle = NA

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
	
    class(fixed_best_angle) <- "numeric";
    if (exists('time_series_dir')) {
        temp = gather_FA_orientation_data(time_series_dir,fixed_best_angle = fixed_best_angle, 
            diagnostic.figure=T);

        temp = gather_FA_orientation_data(time_series_dir,fixed_best_angle = fixed_best_angle,
            min.ratio=2,output_file='FA_orientation_ratio2.Rdata');
        temp = gather_FA_orientation_data(time_series_dir,fixed_best_angle = fixed_best_angle,
            min.ratio=4,output_file='FA_orientation_ratio4.Rdata');
        temp = gather_FA_orientation_data(time_series_dir,fixed_best_angle = fixed_best_angle,
            min.ratio=5,output_file='FA_orientation_ratio5.Rdata');
        temp = gather_FA_orientation_data(time_series_dir,fixed_best_angle = fixed_best_angle,
            min.ratio=6,output_file='FA_orientation_ratio6.Rdata');
    }
}
