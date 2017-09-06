function find_focal_adhesions_full_exp(exp_folder,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overall_start = tic;

i_p = inputParser;

i_p.addRequired('exp_folder',@(x)exist(x,'dir') == 7);

%Adhesion filtering parameters
i_p.addParamValue('i_nums',@(x)isnumeric(x));

i_p.addParamValue('max_adhesion_size',Inf,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('min_adhesion_size',1,@(x)isnumeric(x));
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('min_independent_size',14,@(x)isnumeric(x) && x > 0);

i_p.addParamValue('min_seed_size',NaN,@(x)isnumeric(x) && x > 0);

i_p.addParamValue('no_ad_splitting', 0, @(x) islogical(x) || x == 1 || x == 0);

i_p.addParamValue('max_adhesion_count', Inf, @(x) isnumeric(x));
i_p.addParamValue('stdev_thresh',2,@(x)isnumeric(x));
i_p.addParamValue('per_image_thresh',0,@(x)islogical(x) || x == 0 || x == 1);

i_p.addParamValue('proximity_filter',0,@(x)isnumeric(x) && all(x >= 0));

i_p.addParamValue('confocal_mode',0,@(x)isnumeric(x) && all(x >= 0));

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('paper_figures',0,@(x)x == 1 || x == 0);
i_p.addParamValue('status_messages',1,@(x)x == 1 || x == 0);

i_p.addParamValue('atrous_segmentation',0,@(x)isnumeric(x) && x > 0);
i_p.addParameter('structure_element_size', 5,@(x)isnumeric(x) && x > 1);
i_p.addParameter('atrous_export_level',4,@(x)isnumeric(x) && x > 0);


i_p.parse(exp_folder,varargin{:});

%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');
addpath('../visualize_cell_features');

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%throw out all the parameters that aren't different from the default
%values, we need to do this because the find_focal_adhesion program does
%certain things when variables are not specified, but if we pass the raw
%i_p.Results parameter, it appears that all the parameters are specified,
%even when they were probably mostly default
clean_opts = i_p.Results;
clean_opts = rmfield(clean_opts,'exp_folder');
opt_field_names = fieldnames(clean_opts);
for i = 1:length(opt_field_names)
    if (any(strcmp(opt_field_names{i},i_p.UsingDefaults)))
        clean_opts = rmfield(clean_opts,opt_field_names{i});
    end
end

image_folders = dir(fullfile(exp_folder,'individual_pictures'));
image_folders = image_folders(3:end);
if (not(any(strcmp('i_nums',i_p.UsingDefaults))))
    image_folders = image_folders(i_p.Results.i_nums);
    clean_opts = rmfield(clean_opts,'i_nums');
end

for i = 1:length(image_folders)
    I_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.focal_image);
    try
        if (i_p.Results.atrous_segmentation)
            find_focal_adhesions_atrous(I_file,clean_opts);
        else
            find_focal_adhesions(I_file,clean_opts);
        end
    catch ME
        warning(ME.message);
    end
    disp(['Done with ',I_file]);
end
toc(overall_start);
