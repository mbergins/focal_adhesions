function make_orientation_vis(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('../find_cell_features'));

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

major_axis = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','MajorAxisLength.csv'));
minor_axis = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','MinorAxisLength.csv'));

axis_ratio = major_axis ./ minor_axis;

orientation = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','Orientation.csv'));
dom_angles = csvread(fullfile(exp_dir,'adhesion_props','per_image_dom_angle.csv'));

adj_orientation = abs(adjust_orientation(orientation, dom_angles));

bin_num = 25;
orientation_bins = floor(adj_orientation/(90/bin_num));

orientation_cmap = jet(bin_num);
orientation_cmap = orientation_cmap(bin_num:-1:1,:);

tracking_mat = csvread(fullfile(exp_dir,'tracking_matrices','tracking_seq.csv')) + 1;

output_dir = fullfile(exp_dir,'visualizations','orientation_angle_highlight_jet');
if (not(exist(output_dir,'dir')))
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the Eccen Vis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_num = 1:size(image_dirs)
    this_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
    this_ratio = axis_ratio(:,i_num);
    this_tracking = tracking_mat(:,i_num);
    this_or_bin = orientation_bins(:,i_num);
    
    ad_rows = this_tracking > 0;
    high_eccen_rows = this_ratio >= 3;
    
    high_eccen_rows = ad_rows & high_eccen_rows;
    high_eccen_ads = this_tracking(high_eccen_rows);
    high_eccen_or_bins = this_or_bin(high_eccen_rows);
    or_bin_label = zeros(size(this_data.focal_norm,1),size(this_data.focal_norm,2));
    for i = 1:length(high_eccen_ads)
        ad_num = high_eccen_ads(i);
        or_bin_label(this_data.adhesions_perim == ad_num) = high_eccen_or_bins(i);
    end
    
    low_eccen_rows = ad_rows & not(high_eccen_rows);
    low_eccen_ads = this_tracking(low_eccen_rows);
    low_eccen = ismember(this_data.adhesions,low_eccen_ads);
    
    highlight = create_highlighted_image(this_data.focal_norm,low_eccen,'color_map',[0,0,0]);
    highlight = create_highlighted_image(highlight,or_bin_label,'color_map',orientation_cmap);
    highlight = add_angle_bar(highlight,dom_angles(i_num));
    
    padded_i_num = sprintf('%04d',i_num);
    
    imwrite(highlight,fullfile(output_dir,[padded_i_num,'.png']));
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cor_orientation = adjust_orientation(orientation, dom_angles)

cor_orientation = orientation;

for image_num=1:size(orientation,2)
    this_col = orientation(:,image_num);
    cor_col = apply_new_orientation(this_col,dom_angles(image_num));
    cor_orientation(:,image_num) = cor_col;
end

function new_orientation = apply_new_orientation(orientation,angle)

new_orientation = orientation - angle;

less_neg_ninety = new_orientation < -90;

new_orientation(less_neg_ninety) = new_orientation(less_neg_ninety) + 180;

function highlight = add_angle_bar(highlight,angle,bar_length,edge_spacing)

if(~exist('bar_length','var'))
    bar_length=50;
end

if(~exist('edge_spacing','var'))
    edge_spacing=2;
end

image_size = size(highlight);
mask = lineMask(size(highlight),[(image_size(1)-2),edge_spacing+bar_length], ...
    [image_size(1)-2-bar_length*sind(angle),edge_spacing+bar_length+bar_length*cosd(angle)], ...
    2);

highlight(:,:,1) = highlight(:,:,1) + mask;
highlight(:,:,2) = highlight(:,:,2) + mask;
highlight(:,:,3) = highlight(:,:,3) + mask;
