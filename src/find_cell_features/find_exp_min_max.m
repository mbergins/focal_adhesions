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

min_max = [Inf,-Inf];

temp_image = imread(fullfile(base_dir,image_dirs(1).name,filenames.focal_image));
all_images = zeros(size(temp_image,1),size(temp_image,2),size(image_dirs,1));

for i_num = 1:size(image_dirs,1)
    puncta_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image)));
    if (min(puncta_image(:)) < min_max(1)), min_max(1) = min(puncta_image(:)); end
    if (min(puncta_image(:)) > min_max(2)), min_max(2) = max(puncta_image(:)); end
    
    all_images(:,:,i_num) = puncta_image;
end

output_file = fullfile(base_dir,image_dirs(1).name,filenames.focal_image_min_max);
[output_folder,temp,temp_ext] = fileparts(output_file); %#ok<NASGU>
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end

csvwrite(output_file,min_max);

focal_image_threshold = (2*std(all_images(:)))/range(all_images(:));

csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.focal_image_threshold),focal_image_threshold);

toc;