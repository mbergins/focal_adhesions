function find_cell_mask_full_exp(exp_folder,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overall_start = tic;

i_p = inputParser;

i_p.addRequired('exp_folder',@(x)exist(x,'dir') == 7);

i_p.addParamValue('mask_threshold',0,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('median_filter',0,@(x)x==1 || x==0);
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

for i = 1:length(image_folders)
    mask_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.raw_mask);
    find_cell_mask(mask_file,clean_opts);
    disp(['Done with ',mask_file]);
end
toc(overall_start);
