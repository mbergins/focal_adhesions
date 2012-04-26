function register_with_nifty(exp_dir,varargin)
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
% Collect all the images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%registration will go in a step-wise fashion, register image 2 to image 1,
%then image 3 to 2, then image 4 to 3, ...
for i_num = 1:(size(image_dirs,1)-1)
    FA_file1 = fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image);
    FA_file2 = fullfile(base_dir,image_dirs(i_num+1).name,filenames.focal_image);

    FA_img1 = imread(FA_file1);
    FA_img2 = imread(FA_file2);
        
    img1_nii = make_nii(FA_img1);
    img2_nii = make_nii(FA_img2);
    
    file_nii_temp_1 = fullfile(base_dir,image_dirs(i_num).name,[filenames.focal_image,'.nii']);
    file_nii_temp_2 = fullfile(base_dir,image_dirs(i_num+1).name,[filenames.focal_image,'.nii']);
    
    save_nii(img1_nii,file_nii_temp_1);
    save_nii(img2_nii,file_nii_temp_2);
    
    reg_result_file = fullfile(base_dir,image_dirs(i_num+1).name,'reg_result.nii');
    
    system(['reg_aladin -source ', file_nii_temp_2, ' -target ',file_nii_temp_1, ...
        ' -result ',reg_result_file, ' -aff /dev/null -noRot -rigOnly >/dev/null 2>/dev/null']);
    
    reg_result = load_nii(reg_result_file);
    
    imwrite(reg_result.img,FA_file2);
    
    delete(file_nii_temp_1,file_nii_temp_2,reg_result_file);
    1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
