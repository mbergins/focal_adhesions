function build_resolution_data(varargin)
% BUILD_SIMULATED_DATA    Builds simulated focal adhesion data with known
%                         properties for use in analyzing the quality of
%                         the designed algorithms for studying focal
%                         adhesions
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'BUILD_RESOLUTION_DATA';

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('output_dir', '.', @ischar);

i_p.parse(varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ad_mean_intensity = 0.4237;

background_mean_intensity = 0.1220;

background_noise_var = 0.004;

max_ad_size = 15;
min_ad_size = 1;
ad_padding = ceil(max_ad_size*0.8);
max_ad_separation = 3;
ad_eccentricity = 1;

output_dir = fullfile('..','..','data','simulation','resolution','Images','Paxillin');

image_number = 27;

%Clear out the data directory
png_files = dir(fullfile(output_dir, '*.png'));
for i = 1:size(png_files,1), delete(fullfile(output_dir, png_files(i).name)); end

image_size = [];
for this_i_num = 2:(image_number-1)
    image = [];
    for this_ad_sep = 0:max_ad_separation
        this_row = [];
        
        max_set_size = [2*ad_padding+max_ad_size*ad_eccentricity,2*ad_padding+2*max_ad_size+max_ad_separation];
        for i = min_ad_size:max_ad_size
            temp_ad = fspecial('gaussian',[i*ad_eccentricity,i],(i*ad_eccentricity)/3);
            temp_ad = temp_ad * (1/mean(mean(temp_ad))) * ad_mean_intensity;
            
            temp_ad_set = imnoise(zeros(size(temp_ad,1),2*size(temp_ad,2)+this_ad_sep), ... 
                                  'gaussian',background_mean_intensity,background_noise_var);

            temp_ad_set(1:end,1:size(temp_ad,2)) = temp_ad;
            temp_ad_set(1:end,(size(temp_ad_set,2)-size(temp_ad,2)+1):size(temp_ad_set,2)) = temp_ad;
            
            padded_ad_set = imnoise(zeros(max_set_size),'gaussian',background_mean_intensity,background_noise_var);
            
            row_range = 1:size(temp_ad_set,1);
            row_shift = floor((size(padded_ad_set,1) - size(temp_ad_set,1))/2);
            row_range = row_range + row_shift;

            assert(size(row_range,2) == size(temp_ad_set,1), '%d, %d',size(row_range,2), size(temp_ad_set,1));

            col_range = 1:size(temp_ad_set,2);
            col_shift = floor((size(padded_ad_set,2) - size(temp_ad_set,2))/2);
            col_range = col_range + col_shift;
            assert(size(col_range,2) == size(temp_ad_set,2), '%d, %d',size(col_range,2), size(temp_ad_set,2));

            padded_ad_set(row_range,col_range) = temp_ad_set;
            this_row = [this_row; padded_ad_set]; %#ok<AGROW>
        end
        image = [image,this_row]; %#ok<AGROW>
    end
    image_size = size(image);
    
    if (not(exist('side_image','var')))
        %Create an image to be placed on the side apart from all the other images
        side_ad_size = 3;
        assert(side_ad_size <= max_ad_size, 'Side adhesion size must be less than or equal to the max adhesion size');
        temp_ad = fspecial('gaussian',side_ad_size,side_ad_size/3);
        temp_ad = temp_ad * (ad_mean_intensity/mean(mean(temp_ad)));

        side_image = zeros(size(image,1),2*ad_padding+size(temp_ad,2));
        row_range = 1:size(temp_ad,1);
        row_shift = floor((size(side_image,1) - size(temp_ad,1))/2);
        row_range = row_range + row_shift;

        col_range = 1:size(temp_ad,2);
        col_shift = floor((size(side_image,2) - size(temp_ad,2))/2);
        col_range = col_range + col_shift;

        side_image(row_range,col_range) = temp_ad;
    end
    
    sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
    if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
    imwrite([image, side_image], fullfile(output_dir,[sprintf(sprintf_format,this_i_num), '.png']))
end

%Build the initial blank image
sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
imwrite([zeros(image_size), side_image], fullfile(output_dir,[sprintf(sprintf_format,1), '.png']))

%Build the final blank image
sprintf_format = ['%0', num2str(length(num2str(image_number))), 'd'];
if (not(exist(output_dir,'dir'))); mkdir(output_dir); end
imwrite([zeros(image_size), side_image], fullfile(output_dir,[sprintf(sprintf_format,image_number), '.png']))