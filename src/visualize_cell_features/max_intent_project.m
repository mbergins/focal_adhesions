function max_intent_project(exp_dir,varargin)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the maximum intensity projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_data = read_in_file_set(fullfile(base_dir,image_dirs(1).name),filenames);

all_images = zeros(size(exp_data.focal_image,1),size(exp_data.focal_image,2),size(image_dirs,1));

for i_num = 1:size(image_dirs)
    this_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
    all_images(:,:,i_num) = this_data.focal_image;
end

max_image_projection = double(max(all_images,[],3));
image_set_min_max = csvread(fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image_min_max));
max_image_projection = (max_image_projection - image_set_min_max(1))/range(image_set_min_max);
orig_size = size(max_image_projection);

%Resize the image, to ensure that the adhesion labels can be seen
max_image_projection = imresize(max_image_projection, [1000, NaN]);
new_size = size(max_image_projection);
scale_factor = new_size(1)/orig_size(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather the adhesion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wrap the following reads in trys because the files will sometimes be
%missing, if they are, there weren't any significant events of that class,
%so define placeholder variables
try
    assembly_rows = csvread(fullfile(base_dir,image_dirs(1).name,filenames.assembly_rows));
catch %#ok<CTCH>
    assembly_rows = zeros(0,2);
end

try
    disassembly_rows = csvread(fullfile(base_dir,image_dirs(1).name,filenames.disassembly_rows));
catch %#ok<CTCH>
    disassembly_rows = zeros(0,2);
end

ad_data_set = struct();
ad_data_set.centroid_x = csvread(fullfile(base_dir,image_dirs(1).name,filenames.centroid_x))*scale_factor;
ad_data_set.centroid_y = csvread(fullfile(base_dir,image_dirs(1).name,filenames.centroid_y))*scale_factor;

both_rows = intersect(assembly_rows(:,1),disassembly_rows(:,1));
only_assembly_rows = setdiff(assembly_rows(:,1),both_rows);
only_disassembly_rows = setdiff(disassembly_rows(:,1),both_rows);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlay adhesions on projection and output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_image_projection = add_adhesion_marks(max_image_projection,ad_data_set,both_rows,[1,1,0]);
max_image_projection = add_adhesion_marks(max_image_projection,ad_data_set,only_assembly_rows,[0,1,0]);
max_image_projection = add_adhesion_marks(max_image_projection,ad_data_set,only_disassembly_rows,[1,0,0]);

imwrite(max_image_projection,fullfile(base_dir,image_dirs(1).name,filenames.max_intensity));

all_rows = union(assembly_rows(:,1),disassembly_rows(:,1));

image_size = size(max_image_projection);

add_adhesion_number(fullfile(base_dir,image_dirs(1).name,filenames.max_intensity), ...
    ad_data_set,all_rows,image_size);
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function image = add_adhesion_marks(image,ad_data_set,row_nums,color)

dilation_element = strel('square',5);

for i=1:length(row_nums)
    ad_num = row_nums(i);
    x_pos = round(mean(ad_data_set.centroid_x(ad_num,not(isnan(ad_data_set.centroid_x(ad_num,:))))));
    y_pos = round(mean(ad_data_set.centroid_y(ad_num,not(isnan(ad_data_set.centroid_y(ad_num,:))))));
    
    binary = zeros(size(image));
    binary(x_pos,y_pos) = 1;
    binary = imdilate(binary,dilation_element);
    
    image = create_highlighted_image(image,binary,'color_map',color,'mix_percent',0.65);
end

function add_adhesion_number(image_file,ad_data_set,row_nums,image_size)

for i=1:length(row_nums)
    ad_num = row_nums(i);
    x_pos = round(mean(ad_data_set.centroid_x(ad_num,not(isnan(ad_data_set.centroid_x(ad_num,:))))));
    y_pos = round(mean(ad_data_set.centroid_y(ad_num,not(isnan(ad_data_set.centroid_y(ad_num,:))))));
    
    if (x_pos > image_size(1)*0.95)
        x_pos = image_size(1)*0.95;
    end
    
    if (y_pos > image_size(2)*0.95)
        y_pos = image_size(2)*0.95;
    end
    
    pos_str = [' +',num2str(y_pos),'+',num2str(x_pos)];
    label_str = [' "',num2str(ad_num),'" '];
    command_str = ['convert ', image_file, ' -font VeraBd.ttf -fill ''rgba(0,0,0,0.75)'' -annotate', ...
        pos_str, label_str, ' ', image_file];
    
    system(command_str);
end
