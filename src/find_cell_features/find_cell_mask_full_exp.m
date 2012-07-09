function find_cell_mask_full_exp(exp_folder,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overall_start = tic;

i_p = inputparser;

i_p.addrequired('exp_folder',@(x)exist(x,'dir') == 7);

i_p.addparamvalue('mask_threshold',0,@(x)isnumeric(x) && x > 0);
i_p.addparamvalue('median_filter',0,@(x)x==1 || x==0);
i_p.addparamvalue('single_threshold',0,@(x)x==1 || x==0);
i_p.parse(exp_folder,varargin{:});

%add the folder with all the scripts used in this master program
addpath('matlab_scripts');
addpath('../visualize_cell_features');

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%throw out all the parameters that aren't different from the default
%values, we need to do this because the find_focal_adhesion program does
%certain things when variables are not specified, but if we pass the raw
%i_p.results parameter, it appears that all the parameters are specified,
%even when they were probably mostly default
clean_opts = i_p.results;
clean_opts = rmfield(clean_opts,'exp_folder');
opt_field_names = fieldnames(clean_opts);
for i = 1:length(opt_field_names)
    if (any(strcmp(opt_field_names{i},i_p.usingdefaults)))
        clean_opts = rmfield(clean_opts,opt_field_names{i});
    end
end

image_folders = dir(fullfile(exp_folder,'individual_pictures'));
image_folders = image_folders(3:end);

mask_thresholds = zeros(length(image_folders),1);
for i = 1:length(image_folders)
    mask_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.raw_mask);
    mask_thresholds(i) = find_cell_mask(mask_file,clean_opts);
    disp(['done with ',mask_file]);
end

mask_plot = plot(mask_thresholds);
ylimits = ylim;
ylim([0,ylimits(2)]);
saveas(mask_plot,fullfile(exp_folder,'adhesion_props','cell_mask_thresholds.png'));

if (i_p.results.single_threshold && i_p.Results.mask_threshold > 0)
    clean_opts.mask_threshold = median(mask_thresholds);
    for i = 1:length(image_folders)
        mask_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.raw_mask);
        mask_thresholds(i) = find_cell_mask(mask_file,clean_opts);
    end
end

toc(overall_start);