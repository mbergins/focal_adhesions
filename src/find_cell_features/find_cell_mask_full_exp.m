function find_cell_mask_full_exp(exp_folder,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overall_start = tic;

i_p = inputparser;

i_p.addRequired('exp_folder',@(x)exist(x,'dir') == 7);

i_p.addParamValue('mask_threshold',0,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('median_filter',0,@(x)x==1 || x==0);
i_p.addParamValue('single_threshold',0,@(x)x==1 || x==0);
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
if (not(isempty(strcmp('single_threshold',fieldnames(clean_opts)))))
    clean_opts = rmfield(clean_opts,'single_threshold');
end

image_folders = dir(fullfile(exp_folder,'individual_pictures'));
image_folders = image_folders(3:end);

mask_pix = [];

mask_thresholds = zeros(length(image_folders),1);
for i = 1:length(image_folders)
    mask_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.raw_mask);
    temp = imread(mask_file);
    mask_pix = [mask_pix; temp(:)];
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

single_threshold_val = find_middle_valley(mask_pix);

hist(double(mask_pix(:)),1000);
hold on;
ylimits = ylim;
plot([median(mask_thresholds),median(mask_thresholds)],ylimits,'r');
plot([single_threshold_val,single_threshold_val],ylimits,'g');
print('-dpng', fullfile(exp_folder,'adhesion_props','cell_mask_intensities.png'));
hold off;

mask_plot = plot(mask_thresholds);
ylimits = ylim;
ylim([0,ylimits(2)]);
hold on;
plot([0,length(mask_thresholds)],[median(mask_thresholds),median(mask_thresholds)],'r')
saveas(mask_plot,fullfile(exp_folder,'adhesion_props','cell_mask_thresholds.png'));

if (i_p.Results.single_threshold)
    repro_start = tic;
    disp('Re-processing cell masks with a single threshold');
    clean_opts.mask_threshold = median(mask_thresholds);
    for i = 1:length(image_folders)
        mask_file = fullfile(exp_folder,'individual_pictures',image_folders(i).name,filenames.raw_mask);
        mask_thresholds(i) = find_cell_mask(mask_file,clean_opts);
    end
    toc(repro_start);
end

toc(overall_start);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshold = find_middle_valley(pix_values)

    [heights, intensity] = hist(double(pix_values),1000);
    
    smoothed_heights = smooth(heights,0.01,'loess');
    [~,imax,~,imin]= extrema(smoothed_heights);
    
    %keep in mind that the zmax is sorted by value, so the highest peak is
    %first and the corresponding index is also first in imax, the same pattern
    %hold for zmin and imin
    
    sorted_max_indexes = sort(imax);
    first_max_index = find(sorted_max_indexes == imax(1));
    
    %locate the index between the first two maximums
    min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
    assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));
    threshold = intensity(imin(min_index));