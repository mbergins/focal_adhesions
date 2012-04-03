function find_background_intensity(exp_dir,varargin)
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

upperleft = [];
lowerleft = [];
upperright = [];
lowerright = [];
for i_num = 1:size(image_dirs,1)
    focal_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image)));
    
    corner = focal_image(1:40,1:40);
    upperleft = [upperleft; corner(:)];
    corner = focal_image((size(focal_image,1)-39):end,1:40);
    lowerleft = [lowerleft; corner(:)];
    corner = focal_image(1:40,(size(focal_image,2)-39):end);
    upperright = [upperright; corner(:)];
    corner = focal_image((size(focal_image,1)-39):end,(size(focal_image,2)-39):end);
    lowerright = [lowerright; corner(:)];
end

corner_sets = [upperleft,lowerleft,upperright,lowerright];

min_std_entry = find(std(corner_sets) == min(std(corner_sets)),1,'first');

background_mean = round(mean(corner_sets(:,min_std_entry)));

output_file = fullfile(base_dir,image_dirs(1).name,filenames.background_intensity);
[output_folder,~,~] = fileparts(output_file);
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end
csvwrite(output_file,background_mean);

for i_num = 1:size(image_dirs,1)
    focal_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image)));
    
    corrected_focal_image = focal_image - background_mean;
    corrected_focal_image(corrected_focal_image < 0) = 0;
    
    imwrite(uint16(corrected_focal_image),fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image));    
end

toc;