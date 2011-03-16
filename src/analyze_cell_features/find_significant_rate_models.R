source('FA_analysis_lib.R')
source('bilinear_modeling.R')

args = commandArgs(TRUE);
if (length(args) != 0) {
    debug = FALSE;
    
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
	
    if (exists('exp_dir')) {
        #######################################################################
        # Reading in data files
        #######################################################################
        exp_dirs <- Sys.glob(file.path(exp_dir))
        exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
        
        raw_data = list()
        raw_data$intensity = load_results(exp_dirs,file.path('Average_adhesion_signal.Rdata'));
        
        #######################################################################
        # Filtering and Output
        #######################################################################
        only_signif = filter_results(raw_data$intensity, min.r.sq = -Inf, max.p.val = 0.05)
        
        output_phase_lengths_from_filtered(only_signif,raw_data)
    } else {
        print(paste('The directory that contains all the experiments of interest', 
            'must be specified on the R command line in parameter exp_dir'))
    }
}

