###############################################################################
# FA Orientation Models
###############################################################################

gather_FA_orientation_data <- function(exp_dir, fixed_best_angle = NA) {
    data_set = read_in_orientation_data(exp_dir, min_eccen=3);
    data_set$angle_search = test_dom_angles(data_set$subseted_data$orientation);

    min_angle_indexes = which(min(data_set$angle_search$angle_sds) == data_set$angle_search$angle_sds)

    if (is.na(fixed_best_angle)) {
        data_set$best_angle = median(data_set$angle_search$test_angles[min_angle_indexes])
    } else {
        data_set$best_angle = fixed_best_angle
        data_set$actual_best_angle = median(data_set$angle_search$test_angles[min_angle_indexes])
    }
    data_set$corrected_orientation = apply_new_orientation(data_set$subseted_data$orientation,
        data_set$best_angle)
    
    per_image_dom_angle = find_per_image_dom_angle(data_set$mat, min_eccen=3)

    #diagnostic figure
    pdf(file.path(exp_dir,'..','adhesion_orientation.pdf'))
    layout(rbind(c(1,2),c(3,4)))
    par(bty='n', mar=c(4,4.2,2,0),mgp=c(2,1,0))
    
    hist(data_set$subseted_data$orientation,main='Pos X-axis Reference',
        xlab=paste('Angle n=',dim(data_set$subseted_data)[1],sep=''));
    
    plot(data_set$angle_search$test_angles,data_set$angle_search$angle_sds,typ='l',
        xlab='Dominant Search Angle',ylab='STDEV FA Orientation');
    if (! is.na(fixed_best_angle)) {
        segments(fixed_best_angle,0,fixed_best_angle,2000,col='red')
        segments(data_set$actual_best_angle,0,data_set$actual_best_angle,2000,col='green')
    }

    hist(data_set$corrected_orientation,
        main=paste('Rotated ',data_set$best_angle,'\u00B0 / ',
            sprintf('%0.1f',90-sd(data_set$corrected_orientation)),' STDEV',sep=''),
        xlab=paste('Angle n=',dim(data_set$subseted_data)[1],sep=''));
    
    plot(per_image_dom_angle,xlab='Image Number',ylab='Dominant Angle',ylim=c(0,180));
    graphics.off()

    save(data_set,file=file.path(exp_dir,'..','FA_orientation.Rdata'))
    return(data_set)
}

read_in_orientation_data <- function(time_series_dir,min_eccen = 3) {
    data_set = list();
    data_set$mat$orientation = read.csv(file.path(time_series_dir,'Orientation.csv'),header=F);
    data_set$mat$area = read.csv(file.path(time_series_dir,'Area.csv'),header=F);
    major_axis = read.csv(file.path(time_series_dir,'MajorAxisLength.csv'),header=F);
    minor_axis = read.csv(file.path(time_series_dir,'MinorAxisLength.csv'),header=F);
    
    data_set$mat$eccentricity = major_axis/minor_axis;
    
    

    unlist_data_set = list()
    for (i in names(data_set$mat)) {
        unlist_data_set[[i]] = unlist(data_set$mat[[i]]);
    }
    unlist_data_set = as.data.frame(unlist_data_set);
    data_set$subseted_data = subset(unlist_data_set, !is.nan(unlist_data_set$eccentricity) & 
        unlist_data_set$eccentricity >= min_eccen);
    
    return(data_set);
}

find_per_image_dom_angle <- function(mat_data, min_eccen=3) {
    best_angles = c()
    for (i_num in 1:dim(mat_data$orientation)[2]) {
        this_orientation = mat_data$orientation[,i_num];
        this_eccen = mat_data$eccentricity[,i_num];

        good_ads = !is.na(this_orientation) & this_eccen >= min_eccen;
        
        angle_search = test_dom_angles(this_orientation[good_ads]);
        
        this_best = median(angle_search$test_angles[which(min(angle_search$angle_sds) == angle_search$angle_sds)])
        best_angles = c(best_angles, this_best);
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

    new_angle_sd = c()
    for (angle in angles_to_test) {
        new_orientation = apply_new_orientation(orientation,angle);
        new_angle_sd = c(new_angle_sd, 90-sd(new_orientation));
    }
        
    results$y = new_angle_sd;
    results$angle_sds = new_angle_sd;

    return(results)
}

apply_new_orientation <- function(orientation_data,angle) {
    orientation_data = orientation_data - angle;
    
    less_neg_ninety = orientation_data < -90;
    orientation_data[less_neg_ninety] = orientation_data[less_neg_ninety] + 180;

    return(orientation_data)
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
        temp = gather_FA_orientation_data(time_series_dir,fixed_best_angle = fixed_best_angle);
    }
}
