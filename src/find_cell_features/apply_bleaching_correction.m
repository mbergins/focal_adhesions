function apply_bleaching_correction(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);
i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine single image folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_dir = fullfile(exp_dir, 'individual_pictures');

single_image_folders = dir(image_dir);

assert(strcmp(single_image_folders(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(single_image_folders(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(single_image_folders(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

single_image_folders = single_image_folders(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the Focal Image Bleaching Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

correct_bleaching(image_dir,single_image_folders,filenames.focal_image,exp_dir);
correct_bleaching(image_dir,single_image_folders,filenames.kinase,exp_dir);

toc;

end

function correct_bleaching(image_dir, single_image_folders,file_to_correct,exp_dir)

%Test image and uncorrected files
base_image_file = fullfile(image_dir,single_image_folders(1).name,file_to_correct);
uncorrected_image_file = fullfile(image_dir,single_image_folders(1).name,['uncorrected_',file_to_correct]);

if (not(exist(base_image_file,'file')))
    disp(['Couldn''t find the image file: ',base_image_file]);
    return;
end

if (exist(uncorrected_image_file,'file'))
    disp(['Found the uncorrected image file, skipping image correction:',uncorrected_image_file]);
    return;
end

expression_levels = zeros(1,length(single_image_folders));
after_correction = zeros(1,length(single_image_folders));
sum_expression_levels = zeros(1,length(single_image_folders));
for i=1:length(single_image_folders)
    base_image_file = fullfile(image_dir,single_image_folders(i).name,file_to_correct);
    uncorrected_image_file = fullfile(image_dir,single_image_folders(i).name,['uncorrected_',file_to_correct]);
    if (exist(uncorrected_image_file,'file'))
        disp('Found corrected image file, not applying correction.');
        continue;
    end
    
    image = imread(base_image_file);
    
    expression_levels(i) = mean(image(:));
    sum_expression_levels(i) = sum(image(:));
    
    movefile(base_image_file,uncorrected_image_file);
    
    cor_image = image.*(expression_levels(1)/expression_levels(i));
    after_correction(i) = mean(cor_image(:));
    imwrite(cor_image, base_image_file,'Bitdepth',16);
end

summary_stats_dir = fullfile(exp_dir,'adhesion_props','photo_bleaching_correction');
if (not(exist(summary_stats_dir,'dir')))
    mkdir(summary_stats_dir);
end

csvwrite(fullfile(summary_stats_dir,[file_to_correct,'_mean_levels.csv']),expression_levels)
csvwrite(fullfile(summary_stats_dir,[file_to_correct,'_integrated_levels.csv']),sum_expression_levels)

end