function build_thresholded_set(I_file,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.parse(I_file,varargin{:});

%read in and normalize the input focal adhesion image

focal_image  = double(imread(I_file));
if (length(size(focal_image)) > 2) 
    focal_image = focal_image(:,:,1);
end
focal_image_norm = (focal_image - min(focal_image(:)))/range(focal_image(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_filt = fspecial('disk',11);
blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
high_passed_image = focal_image - blurred_image;

mean_val = mean(high_passed_image(:));
stdev_val = std(high_passed_image(:));

stdev_intervals = 1:0.5:4;

% stdev_intervals = linspace(-3,3,20);
% temp = jet(floor(length(stdev_intervals)/2));
% c_map = jet(ceil(length(stdev_intervals)/2));
% c_map = [temp(end:-1:1,:);c_map];
thresholds = mean_val + stdev_val*stdev_intervals;

%collected the thresholded images
threshed_images = cell(0);
threshed_images_perims = cell(0);
% threshed_image_overlay = ones([size(high_passed_image),3]);
% i_pix_count = length(high_passed_image(:));
% c_map = jet(length(stdev_intervals));
for i = 1:length(thresholds)
    threshed_images{i} = find_threshed_image(high_passed_image,thresholds(i));

    %identify and remove adhesions on the immediate edge of the image
%     threshed_images{i} = remove_edge_adhesions(threshed_images{i});
    
    threshed_images_perims{i} = bwperim(threshed_images{i});
    1;
%     threshed_image_overlay(find(threshed_images{i})) = c_map(i,1);
%     threshed_image_overlay(find(threshed_images{i})+i_pix_count) = c_map(i,2);
%     threshed_image_overlay(find(threshed_images{i})+i_pix_count*2) = c_map(i,3);    
end

%create highlighted images
highlighted_images = cell(0);
for i = 1:length(threshed_images)
    highlighted_images{i} = create_highlighted_image(focal_image_norm,threshed_images_perims{i},'color_map',[1,1,0]);
end

%file output
[folder,base_name] = fileparts(I_file);
imwrite(focal_image_norm, fullfile(folder,[base_name,'_norm.png']));
output_files = cell(0);
for i=1:length(highlighted_images)
    output_files{i} = fullfile(folder,[base_name,'_',num2str(stdev_intervals(i)),'.png']);
    imwrite(highlighted_images{i},output_files{i});
end

%image labeling
font_size = num2str(round((30/520)*size(focal_image,1)));
command = ['convert ', fullfile(folder,[base_name,'_norm.png']), ' -font Ubuntu-B.ttf -pointsize ', font_size, ...
    ' -gravity southwest -fill white -annotate 0 "Original Image" ', ... 
    fullfile(folder,[base_name,'_norm.png'])];
system(command);
for i=1:length(highlighted_images)
    command = ['convert ', output_files{i}, ' -font Ubuntu-B.ttf -pointsize ', font_size, ... 
        ' -gravity southwest -fill white -annotate 0 "Stdev Threshold: ', ...
        num2str(stdev_intervals(i)), '" ', output_files{i}];
    system(command);
end

%build the thumb files
thumbs_folder = fullfile(folder,'thumbs');
mkdir(thumbs_folder);
system(['mogrify -format jpg -path ',thumbs_folder,' -thumbnail 400 ',folder,'/*.png']);

toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshed_image = find_threshed_image(high_passed_image, filter_thresh)

if (length(filter_thresh) == 1)
    threshed_image = high_passed_image >= filter_thresh;
else
    high_threshed_image = high_passed_image >= filter_thresh(2);
    high_threshed_image = remove_edge_adhesions(high_threshed_image);
    
    low_threshed_image = high_passed_image >= filter_thresh(1);
    low_thresh_bwlabel = bwlabel(low_threshed_image,4);
    
    overlap_labels = unique(low_thresh_bwlabel.*high_threshed_image);
    if (overlap_labels(1) == 0)
        overlap_labels = overlap_labels(2:end);
    end
    
    threshed_image = ismember(low_thresh_bwlabel,overlap_labels);
end
end

function cleaned_binary = remove_edge_adhesions(threshed_image,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'REMOVE_EDGE_ADHESIONS';

i_p.addRequired('threshed_image',@(x)isnumeric(x) || islogical(x));
i_p.addParamValue('binary_shift',0,@islogical);

i_p.parse(threshed_image,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edge_border = ones(size(threshed_image));
if (not(any(strmatch('binary_shift',i_p.UsingDefaults))))
    edge_border = bwperim(i_p.Results.binary_shift);
else
    edge_border = bwperim(edge_border);
end

[row_indexes,col_indexes] = ind2sub(size(threshed_image), find(edge_border));
edge_adhesions = bwselect(threshed_image,col_indexes,row_indexes,4);

cleaned_binary = threshed_image & not(edge_adhesions);
end

function high_image = create_highlighted_image(I,high,varargin)
%CREATE_HIGHLIGHTED_IMAGE    add highlights to an image
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS) adds green highlights to
%   image 'I', using the binary image 'HIGHLIGHTS' as the guide
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS,'color_map',[R,G,B]) adds
%   highlights of color specified by the RGB sequence '[R,G,B]' to image
%   'I', using the binary image 'HIGHLIGHTS' as the guide

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'CREATE_HIGHLIGHTED_IMAGE';

i_p.addRequired('I',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('high',@(x)(isnumeric(x) || islogical(x)));

i_p.parse(I,high);

i_p.addParamValue('color_map',jet(double(max(high(:)))),@(x)(all(high(:) == 0) || (isnumeric(x) && (size(x,1) >= max(unique(high))))));
i_p.addParamValue('mix_percent',1,@(x)(isnumeric(x)));

i_p.parse(I,high,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_size = size(I);

if (size(image_size) < 3)
    high_image_red = I;
    high_image_green = I;
    high_image_blue = I;
else
    high_image_red = I(:,:,1);
    high_image_green = I(:,:,2);
    high_image_blue = I(:,:,3);
end

if (all(high(:) == 0))
    high_image = cat(3,high_image_red,high_image_green,high_image_blue);
    return
end

labels = unique(high);
assert(labels(1) == 0)
labels = labels(2:end);

for i=1:length(labels)
    indexes = high == labels(i);
    
    this_cmap = i_p.Results.color_map(labels(i),:);
    
    high_image_red(indexes) = this_cmap(1)*i_p.Results.mix_percent + high_image_red(indexes)*(1-i_p.Results.mix_percent);
    high_image_green(indexes) = this_cmap(2)*i_p.Results.mix_percent + high_image_green(indexes)*(1-i_p.Results.mix_percent);
    high_image_blue(indexes) = this_cmap(3)*i_p.Results.mix_percent + high_image_blue(indexes)*(1-i_p.Results.mix_percent);
end

high_image = cat(3,high_image_red,high_image_green,high_image_blue);

end

% function ad_perims = find_ad_outlines(adhesions)
% 
% ad_zamir_perim = zeros(size(ad_zamir));
% for i = 1:max(ad_zamir(:))
%     assert(any(any(ad_zamir == i)), 'Error: can''t find ad number %d', i);
%     this_ad = zeros(size(ad_zamir));
%     this_ad(ad_zamir == i) = 1;
%     ad_zamir_perim(bwperim(this_ad)) = i;
% end
% 
