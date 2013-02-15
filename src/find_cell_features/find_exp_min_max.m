function find_exp_min_max(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath('matlab_scripts');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FA image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_images = [];
for i_num = 1:size(image_dirs,1)
    image_file_name = fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image);
    
    if (any(size(all_images) == 0))
        temp_image = double(imread(image_file_name));
        all_images = zeros(size(temp_image,1),size(temp_image,2),size(image_dirs,1));
        all_images(:,:,i_num) = temp_image;
    else
        all_images(:,:,i_num) = double(imread(image_file_name)); %#ok<AGROW>
    end
end

trimmed_vals = trim_data_set(all_images(:),1E-4);
fa_min_max = [trimmed_vals(1),trimmed_vals(end)];

csvwrite_with_folder_creation(fullfile(base_dir,image_dirs(1).name,filenames.focal_image_min_max), ...
    fa_min_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinase Min/Max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_name = fullfile(base_dir,image_dirs(1).name,filenames.kinase);

if (exist(image_file_name,'file') == 2)
    all_images = zeros(all_images);
    for i_num = 1:size(image_dirs,1)
        image_file_name = fullfile(base_dir,image_dirs(i_num).name,filenames.kinase);
        all_images(:,:,i_num) = double(imread(image_file_name));
    end
    
    trimmed_vals = trim_data_set(all_images(:),1E-4);
    kinase_min_max = [trimmed_vals(1),trimmed_vals(end)];
    
    csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.kinase_min_max),...
        kinase_min_max);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell Mask Min/Max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_name = fullfile(base_dir,image_dirs(1).name,filenames.raw_mask);

if (exist(image_file_name,'file') == 2)
    all_images = zeros(size(temp_image,1),size(temp_image,2),size(image_dirs,1));
    for i_num = 1:size(image_dirs,1)
        image_file_name = fullfile(base_dir,image_dirs(i_num).name,filenames.raw_mask);
        all_images(:,:,i_num) = double(imread(image_file_name));
    end
    
    trimmed_vals = trim_data_set(all_images(:),1E-4);
    raw_mask_min_max = [trimmed_vals(1),trimmed_vals(end)];
    
    csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.raw_mask_min_max),...
        raw_mask_min_max);
end

toc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trimmed_set = trim_data_set(data,fraction)

trimmed_set = sort(data);
remove_limit = round(length(trimmed_set)*fraction);

trimmed_set = trimmed_set(remove_limit:(end - remove_limit));

end