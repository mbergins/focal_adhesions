function make_eccen_filtered_vis(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

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

major_axis = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','MajorAxisLength.csv'));
minor_axis = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','MinorAxisLength.csv'));

axis_ratio = major_axis ./ minor_axis;

tracking_mat = csvread(fullfile(exp_dir,'tracking_matrices','tracking_seq.csv')) + 1;

output_dir = fullfile(exp_dir,'visualizations','axis_ratio_highlight');
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
    
    ad_rows = this_tracking > 0;
    high_eccen_rows = this_ratio >= 3;
    
    high_eccen_rows = ad_rows & high_eccen_rows;
    low_eccen_rows = ad_rows & not(high_eccen_rows);
    
    high_eccen_ads = this_tracking(high_eccen_rows);
    low_eccen_ads = this_tracking(low_eccen_rows);
    
    high_eccen_filled = ismember(this_data.adhesions,high_eccen_ads);
    low_eccen_filled = ismember(this_data.adhesions,low_eccen_ads);
    
    highlight = create_highlighted_image(this_data.focal_norm,high_eccen_filled,'color_map',[0,1,0],'mix_percent',1);
    highlight = create_highlighted_image(highlight,low_eccen_filled,'color_map',[1,0,0],'mix_percent',1);
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

x_pos = round(centroid(1));
y_pos = round(centroid(2));

binary = zeros(size(image));
binary(x_pos,y_pos) = 1;
binary = imdilate(binary,dilation_element);

image = create_highlighted_image(image,binary,'color_map',color);
