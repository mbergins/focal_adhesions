function find_exp_thresholds(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('stdev_thresh',4,@(x)isnumeric(x) && x > 0);
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

temp_image = imread(fullfile(base_dir,image_dirs(1).name,filenames.focal_image));
all_images = zeros(size(temp_image,1),size(temp_image,2),size(image_dirs,1));
all_high_passed = zeros(size(temp_image,1),size(temp_image,2),size(image_dirs,1));

for i_num = 1:size(image_dirs,1)
    puncta_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image)));
    
    I_filt = fspecial('disk',11);
    blurred_image = imfilter(puncta_image,I_filt,'same',mean(puncta_image(:)));
    high_passed_image = puncta_image - blurred_image;
    
    all_images(:,:,i_num) = puncta_image;
    all_high_passed(:,:,i_num) = high_passed_image;
    disp([i_num,size(image_dirs,1)])
end
min_max = [min(all_images(:)),max(all_images(:))];

output_file = fullfile(base_dir,image_dirs(1).name,filenames.focal_image_min_max);
[output_folder,temp,temp_ext] = fileparts(output_file); %#ok<NASGU>
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end

csvwrite(output_file,min_max);

% im_std = std(all_images,1,3);
% im_mean = mean(all_images,3);
% im_median = median(all_images,3);
% mean_med_diff = im_median - im_mean;
% ans_mean = mean(all_images_anscomb,3);
% ans_var = var(all_images_anscomb,0,3);
% threshed = mean_med_diff > -2 & mean_med_diff < 2;
% im_cv = im_std./im_mean;

focal_image_threshold = mean(all_high_passed(:)) + i_p.Results.stdev_thresh*std(all_high_passed(:));
csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.focal_image_threshold),focal_image_threshold);

toc;