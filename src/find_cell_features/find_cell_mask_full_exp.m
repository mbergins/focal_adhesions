function find_cell_mask_full_exp(exp_folder,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overall_start = tic;

i_p = inputParser;

i_p.addRequired('exp_folder',@(x)exist(x,'dir') == 7);

i_p.addParamValue('mask_threshold',0,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('median_filter',0,@(x)x==1 || x==0);
i_p.addParamValue('single_threshold',0,@(x)x==1 || x==0);

i_p.addParamValue('debug',0,@(x)x==1 || x==0);

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
clean_opts = i_p.Results;
clean_opts = rmfield(clean_opts,'exp_folder');
opt_field_names = fieldnames(clean_opts);
for i = 1:length(opt_field_names)
    if (any(strcmp(opt_field_names{i},i_p.UsingDefaults)))
        clean_opts = rmfield(clean_opts,opt_field_names{i});
    end
end

%remove single_threshold parameter, if present
if (any(strcmp('single_threshold',fieldnames(clean_opts))))
    clean_opts = rmfield(clean_opts,'single_threshold');
end

image_folders = dir(fullfile(exp_folder,'individual_pictures'));
image_folders = image_folders(3:end);

mask_thresholds = zeros(length(image_folders),1);
for i = 1:length(image_folders)
    mask_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.raw_mask);
    mask_thresholds(i) = find_cell_mask(mask_file,clean_opts);
    if (mod(i,10) == 0)
        disp(['Done with ',mask_file]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (not(exist(fullfile(exp_folder,'adhesion_props'),'dir')))
    mkdir(fullfile(exp_folder,'adhesion_props'))
end

mask_plot = plot(mask_thresholds);
ylimits = ylim;
ylim([0,ylimits(2)]);
hold on;
plot([0,length(mask_thresholds)],[nanmedian(mask_thresholds),nanmedian(mask_thresholds)],'r')
saveas(mask_plot,fullfile(exp_folder,'adhesion_props','image_analysis','cell_mask_thresholds.png'));
hold off;

if (i_p.Results.single_threshold)
    repro_start = tic;
    disp('Re-processing cell masks with a single threshold');
    clean_opts.mask_threshold = nanmedian(mask_thresholds);
    for i = 1:length(image_folders)
        mask_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.raw_mask);
        mask_thresholds(i) = find_cell_mask(mask_file,clean_opts);
    end
    toc(repro_start);
end

toc(overall_start);