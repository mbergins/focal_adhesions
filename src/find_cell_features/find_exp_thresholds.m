function find_exp_thresholds(exp_dir,varargin)
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

temp_image = imread(fullfile(base_dir,image_dirs(1).name,filenames.focal_image));
all_images = zeros(size(temp_image,1),size(temp_image,2),size(image_dirs,1));
all_high_passed = zeros(size(temp_image,1),size(temp_image,2),size(image_dirs,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect all the images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_num = 1:size(image_dirs,1)
    puncta_image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image)));            
    
    I_filt = fspecial('disk',11);
    blurred_image = imfilter(puncta_image,I_filt,'same',mean(puncta_image(:)));
    high_passed_image = puncta_image - blurred_image;
    
    all_images(:,:,i_num) = puncta_image;
    all_high_passed(:,:,i_num) = high_passed_image;
    if (mod(i_num,10) == 0)
        disp([i_num,size(image_dirs,1)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine min/max and identification thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trimmed_pix_values = trim_data_set(all_images(:),1E-4);
min_max = [trimmed_pix_values(1),trimmed_pix_values(end)];

output_file = fullfile(base_dir,image_dirs(1).name,filenames.focal_image_min_max);
[output_folder,~,~] = fileparts(output_file);
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end
csvwrite(output_file,min_max);

high_pass_mean = mean(all_high_passed(:));
high_pass_std = std(all_high_passed(:));
csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.focal_image_threshold),...
    [high_pass_mean,high_pass_std]);

per_image_threshold = zeros(size(image_dirs,1),2);
for i_num = 1:size(image_dirs,1)
    temp = all_high_passed(:,:,i_num);
    temp = trim_data_set(temp(:),5E-6);

    per_image_threshold(i_num,:) = [mean(temp(:)), std(temp(:))];
end

%diagnostic diagram
hist(all_high_passed(:),100);
xlabel('High Pass Filtered Intensity','FontSize',16,'FontName','Helvetica');
ylabel('Pixel Count','FontSize',16,'FontName','Helvetica');
y_limits = ylim();

for i=1:4
    this_thresh = high_pass_mean+high_pass_std*i;
    line([this_thresh,this_thresh],[0,y_limits(2)],'Color','red', ... 
        'LineStyle','--','LineWidth',3);
end

set(gca, 'FontName','Helvetica','FontSize',16,'Box','off');
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', fullfile(base_dir,image_dirs(1).name,filenames.focal_image_threshold_plot));
close;

toc;
end

function trimmed_set = trim_data_set(data,percentage)

trimmed_set = sort(data);
remove_limit = round(length(trimmed_set)*percentage);

trimmed_set = trimmed_set(remove_limit:(end - remove_limit));

end