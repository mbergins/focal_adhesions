function find_exp_thresholds(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('stdev_thresh',4,@(x)isnumeric(x) && all(x > 0));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

stdev_thresh = i_p.Results.stdev_thresh;

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
    if (mod(i_num,10) == 0)
        disp([i_num,size(image_dirs,1)])
    end
end
sorted_pix_vals = sort(all_images(:));
remove_limit = round(length(sorted_pix_vals)/100000);
sorted_pix_vals = sorted_pix_vals(remove_limit:(length(sorted_pix_vals)-remove_limit));
min_max = [min(sorted_pix_vals),max(sorted_pix_vals)];

output_file = fullfile(base_dir,image_dirs(1).name,filenames.focal_image_min_max);
[output_folder,~,~] = fileparts(output_file);
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end

csvwrite(output_file,min_max);

focal_image_threshold = mean(all_high_passed(:)) + i_p.Results.stdev_thresh*std(all_high_passed(:));
csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.focal_image_threshold),focal_image_threshold);

hist(all_high_passed(:),100);
xlabel('High Pass Filtered Intensity','FontSize',16,'FontName','Helvetica');
ylabel('Pixel Count','FontSize',16,'FontName','Helvetica');
y_limits = ylim();
line([focal_image_threshold,focal_image_threshold],[0,y_limits(2)],'Color','red', ... 
    'LineStyle','--','LineWidth',3);
set(gca, 'FontName','Helvetica','FontSize',16,'Box','off');
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', fullfile(base_dir,image_dirs(1).name,filenames.focal_image_threshold_plot));
close;

toc;