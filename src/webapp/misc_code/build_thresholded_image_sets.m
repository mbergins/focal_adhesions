function build_thresholded_image_sets(I_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load image;
tic;

focal_image  = double(imread(I_file));
if (length(size(focal_image)) > 2) 
    focal_image = focal_image(:,:,1);
end

focal_image = (focal_image - min(focal_image(:)))/range(focal_image(:));

focal_image_norm = (focal_image - min(focal_image(:)))/range(focal_image(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_filt = fspecial('disk',11);
blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
high_passed_image = focal_image - blurred_image;

mean_val = mean(high_passed_image(:));
stdev_val = std(high_passed_image(:));

stdev_intervals = 1:0.5:6;

thresholds = mean_val + stdev_val*stdev_intervals;

%collected the thresholded images
threshed_images = cell(0);
threshed_images_perims = cell(0);
for i = 1:length(thresholds)
    threshed_images{i} = find_threshed_image(high_passed_image,thresholds(i));

    threshed_images_perims{i} = bwperim(threshed_images{i});
end

%create highlighted images
highlighted_images = cell(0);
for i = 1:length(threshed_images)
    highlighted_images{i} = create_highlighted_image(focal_image_norm,threshed_images_perims{i},[1,1,0]);
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
command = ['convert ', fullfile(folder,[base_name,'_norm.png']), ' -weight Bold -pointsize ', font_size, ...
    ' -gravity southwest -fill white -annotate 0 "Original Image" ', ... 
    fullfile(folder,[base_name,'_norm.png'])];
system(command);
for i=1:length(highlighted_images)
    command = ['convert ', output_files{i}, ' -weight Bold  -pointsize ', font_size, ... 
        ' -gravity southwest -fill white -annotate 0 "stdev Threshold: ', ...
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

function cleaned_binary = remove_edge_adhesions(threshed_image)

edge_border = ones(size(threshed_image));
edge_border = bwperim(edge_border);

[row_indexes,col_indexes] = ind2sub(size(threshed_image), find(edge_border));
edge_adhesions = bwselect(threshed_image,col_indexes,row_indexes,4);

cleaned_binary = threshed_image & not(edge_adhesions);
end

function high_image = create_highlighted_image(I,high,color_map)

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
    
    this_cmap = color_map(labels(i),:);
    
    high_image_red(indexes) = this_cmap(1);
    high_image_green(indexes) = this_cmap(2);
    high_image_blue(indexes) = this_cmap(3);
end

high_image = cat(3,high_image_red,high_image_green,high_image_blue);

end
