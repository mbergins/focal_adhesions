function make_filtered_vis(exp_dir,varargin)
%  MAKE_FILTERED_VIS    Creates visualization of FA highlighted by various
%                       properties

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('type','ratio',@ischar);
i_p.addParamValue('min_value',1,@isnumeric);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>
image_dirs = image_dirs(3:end);

tracking_mat = csvread(fullfile(exp_dir,'tracking_matrices','tracking_seq.csv')) + 1;

%highlight_decision: 1 for adhesion passed the filter, 0 for didn't pass
%the filter
%   -seed the matrix with all of the positions in the matrix with adhesions
%   present
highlight_decision = zeros(size(tracking_mat));
highlight_decision(tracking_mat > 0) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the type parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch i_p.Results.type
    case 'ratio'
        major_axis = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','MajorAxisLength.csv'));
        minor_axis = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','MinorAxisLength.csv'));
        
        axis_ratio = major_axis ./ minor_axis;
        
        highlight_decision(axis_ratio < 2) = 2;
        highlight_decision(axis_ratio >= 2) = 1;
        
        output_dir = fullfile(exp_dir,'visualizations','axis_ratio_highlight');
    case 'lifetime'
        area = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','Area.csv'));
        lifetime = sum(not(isnan(area)),2);
        
        highlight_decision(lifetime < i_p.Results.min_value,:) = 0;
        
        output_dir = fullfile(exp_dir,'visualizations',['lifetime_',num2str(i_p.Results.min_value)]);
    case 'FA_dist'
        FA_dist = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','Dist_to_FA_cent.csv'));
        FA_dist_mean = nanmean(FA_dist,2);
        
        if(any(strcmp('min_value',i_p.UsingDefaults)))
            p_tile = 0.40;
        else
            p_tile = i_p.Results.min_value;
        end
        
        min_val = quantile(FA_dist_mean,p_tile);
        
        highlight_decision(FA_dist_mean < min_val,:) = 0;
        
        output_dir = fullfile(exp_dir,'visualizations',['FA_dist_',num2str(p_tile)]);
    case 'FA_CHull_dist'
        FA_dist = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','CHull_dist.csv'));
        FA_dist_mean = nanmean(FA_dist,2);
        
        if(any(strcmp('min_value',i_p.UsingDefaults)))
            p_tile = 0.40;
        else
            p_tile = i_p.Results.min_value;
        end
        
        max_val = quantile(FA_dist_mean,p_tile);
        
        highlight_decision(FA_dist_mean > max_val,:) = 0;
        
        output_dir = fullfile(exp_dir,'visualizations',['FA_CHull_dist_',num2str(p_tile)]);
    case 'FA_angle'
        FA_angle = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','FA_angle_recentered.csv'));
        FA_angle_mean = nanmean(abs(FA_angle),2);
        
        max_val = 90;
        
        %mark adhesions above the maximum angle - on a per image basis
        highlight_decision(abs(FA_angle) > max_val) = 2;
        
        %mark adhesions above the maximum angle - on an average angle basis
%         highlight_decision(FA_angle_mean < max_val,:) = 2;
        
        output_dir = fullfile(exp_dir,'visualizations',['FA_angle_',num2str(max_val)]);
    case 'FA_angle_outside'
        FA_angle = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','FA_angle_recentered.csv'));
        FA_angle_mean = nanmean(abs(FA_angle),2);
        
        min_val = 60;
        max_val = 120;
        
        %remove adhesions above the maximum angle
        highlight_decision(FA_angle_mean <= min_val,:) = 0;
        highlight_decision(FA_angle_mean > max_val,:) = 0;
        
        output_dir = fullfile(exp_dir,'visualizations',['FA_angle_outside_',num2str(60)]);
    case 'zaozao'
        FA_angle = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','FA_angle_recentered.csv'));
        FA_angle_mean = nanmean(abs(FA_angle),2);
        
        FA_dist = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','CHull_dist.csv'));
        FA_dist_mean = nanmean(FA_dist,2);
        
        dist_cutoff = quantile(FA_dist_mean,0.4);
        
        highlight_decision(FA_angle_mean <= 80 & FA_dist_mean <= dist_cutoff,:) = 1;
        highlight_decision(FA_angle_mean > 80 & FA_angle_mean <= 120 ...
            & FA_dist_mean <= dist_cutoff,:) = 2;
        highlight_decision(FA_dist_mean > dist_cutoff | FA_angle_mean > 120,:) = 3;
        
        output_dir = fullfile(exp_dir,'visualizations','ZZ');
    case 'live_at_point'
        highlight_decision(tracking_mat(:,15) <= 0,:) = 2;
        highlight_decision(tracking_mat(:,15) > 0,:) = 1;
        
        output_dir = fullfile(exp_dir,'visualizations','live_at_15');
    case 'at_birth'
        highlight_decision(tracking_mat > 0) = 2;
        for ad_num = 1:size(tracking_mat,1) 
            highlight_decision(ad_num,find(tracking_mat(ad_num,:) > 0,1,'first')) = 1;
        end
        output_dir = fullfile(exp_dir,'visualizations','at_birth');        
    otherwise
        disp(['Undefined visualization type requested: "',i_p.Results.type, '" exiting.']);
        return
end

if (not(exist(output_dir,'dir')))
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the Eccen Vis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color_maps = [0,1,0;1,0,0;0,0,1];

for i_num = 1:size(image_dirs)
    this_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
    this_decision = highlight_decision(:,i_num);
    this_tracking = tracking_mat(:,i_num);
    
    ad_rows = this_tracking > 0;
    
    highlight_nums = unique(this_decision(ad_rows));
    
    highlight = this_data.focal_norm;
    for h_num = highlight_nums'
        these_ads = this_tracking(ad_rows & this_decision == h_num);
        
        mask = ismember(this_data.adhesions,these_ads);
        highlight = create_highlighted_image(highlight,mask, ...
            'color_map',color_maps(h_num,:),'mix_percent',1);
    end
    
    if (any(strcmp('adhesion_centroid',fieldnames(this_data))))
        highlight = add_centroid_mark(highlight,this_data.adhesion_centroid,[0,0,1]);
    end
    
    padded_i_num = sprintf('%04d',i_num);
    out_file = fullfile(output_dir,[padded_i_num,'.png']);
    %     highlight = imresize(highlight,[591, NaN]);
    imwrite(highlight,out_file);
    %     system(['convert ',out_file,' -pointsize 36 -gravity southwest -fill "#FFFFFF" -annotate 0 "High FAAI" ',out_file]);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function image = add_centroid_mark(image,centroid,color)

dilation_element = strel('square',5);

x_pos = round(centroid(2));
y_pos = round(centroid(1));

binary = zeros(size(image));
binary(x_pos,y_pos) = 1;
binary = imdilate(binary,dilation_element);

image = create_highlighted_image(image,binary,'color_map',color);
