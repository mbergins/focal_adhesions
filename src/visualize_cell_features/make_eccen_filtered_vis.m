function make_filtered_vis(exp_dir,varargin)

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

eccentricity = major_axis ./ minor_axis;

tracking_mat = csvread(fullfile(exp_dir,'tracking_matrices','tracking_seq.csv')) + 1;

output_dir = fullfile(exp_dir,'visualizations','eccentricity_highlight');
if (not(exist(output_dir,'dir')))
    mkdir(output_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the Eccen Vis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_num = 1:size(image_dirs)
    this_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
    this_eccen = eccentricity(:,i_num);
    this_tracking = tracking_mat(:,i_num);
    
    ad_rows = this_tracking > 0;
    high_eccen_rows = this_eccen > 2;
    
    high_eccen_rows = ad_rows & high_eccen_rows;
    low_eccen_rows = ad_rows & not(high_eccen_rows);
    
    high_eccen_ads = this_tracking(high_eccen_rows);
    low_eccen_ads = this_tracking(low_eccen_rows);
    
    high_eccen_perim = ismember(this_data.adhesions_perim,high_eccen_ads);
    low_eccen_perim = ismember(this_data.adhesions_perim,low_eccen_ads);
    
    highlight = create_highlighted_image(this_data.focal_norm,high_eccen_perim,'color_map',[0,1,0]);
    highlight = create_highlighted_image(highlight,low_eccen_perim,'color_map',[1,0,0]);
    
    padded_i_num = sprintf('%04d',i_num);
    
    imwrite(highlight,fullfile(output_dir,[padded_i_num,'.png']));
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%