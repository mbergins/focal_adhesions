function [varargout] = find_focal_adhesions(I_file,varargin)
% FIND_FOCAL_ADHESIONS    locates the focal adhesions in a given image,
%                         optionally returns the segmented image or writes
%                         the segmented image to a file
%
%   find_focal_adhesions(I,OPTIONS) locate the focal adhesions in image
%   file, 'I'
%
%   Options:
%
%       -cell_mask: file which contains the cell mask
%       -filter_size: size of the averaging filter to use, defaults to 23
%       -output_dir: folder used to hold all the results, defaults to the
%        same folder as the image file, 'I'
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.StructExpand = true;
i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.parse(I_file);

%Adhesion filtering parameters
i_p.addParamValue('min_adhesion_size',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('min_independent_size',14,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('no_ad_splitting', 0, @(x) islogical(x) || x == 1 || x == 0);
i_p.addParamValue('max_adhesion_count', Inf, @(x) isnumeric(x));
i_p.addParamValue('stdev_thresh',2,@(x)isnumeric(x) && all(x > 0));
i_p.addParamValue('proximity_filter',0,@(x)isnumeric(x) && all(x >= 0));

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('paper_figures',0,@(x)x == 1 || x == 0);
i_p.addParamValue('status_messages',1,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');
addpath('../visualize_cell_features');

filenames = add_filenames_to_struct(struct());

%read in the cell mask image if defined in parameter set
if (exist(fullfile(fileparts(I_file),filenames.cell_mask),'file'))
    cell_mask = imread(fullfile(fileparts(I_file),filenames.cell_mask));
end

filter_vals = csvread(fullfile(fileparts(I_file),filenames.focal_image_threshold));
filter_thresh = filter_vals(1) + filter_vals(2)*i_p.Results.stdev_thresh;

%read in and normalize the input focal adhesion image
focal_image  = double(imread(I_file));
image_set_min_max = csvread(fullfile(fileparts(I_file),filenames.focal_image_min_max));
focal_normed = (focal_image - image_set_min_max(1))/(range(image_set_min_max));

output_dir = fileparts(I_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply filter to find adhesion regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_filt = fspecial('disk',i_p.Results.filter_size);
blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
high_passed_image = focal_image - blurred_image;

filter_thresh = [mean(high_passed_image(:)) + std(high_passed_image(:))*i_p.Results.stdev_thresh];

threshed_image = find_threshed_image(high_passed_image,filter_thresh,i_p.Results.proximity_filter);

%identify and remove adhesions on the immediate edge of the image
threshed_image = remove_edge_adhesions(threshed_image);

%filter out small adhesions if requested
if (i_p.Results.min_adhesion_size > 1)
    labeled_thresh = bwlabel(threshed_image,4);
    
    props = regionprops(labeled_thresh,'Area'); %#ok<MRPBW>
    labeled_thresh = ismember(labeled_thresh, find([props.Area] >= i_p.Results.min_adhesion_size));
    
    threshed_image = labeled_thresh > 0;
end

%adding a check for finding adhesions, if didn't find any, output error
%file
if (sum(sum(threshed_image)) == 0)
    no_ad_found_file = fullfile(output_dir, 'no_ads_found.txt');
    system(['touch ', no_ad_found_file]);
    warning('FA:noadsfound','Didn''t find any adhesions');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adhesion Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%If adhesion splitting is off ('no_ad_splitting'), then watershed based
%segmentation isn't needed because all the watershed method adds is the
%ability to split touching adhesions. Also, we need to know the pixel size
%in order to select a threshold for having adhesions remain as seperate
%enties when touching, so we also check for that before using the watershed
%segmentation
if (i_p.Results.no_ad_splitting || any(strcmp('min_independent_size',i_p.UsingDefaults)))
    %if splitting is off, there is no need to use the fancy watershed based
    %segmentation methods, just identify the connected areas
    ad_zamir = bwlabel(threshed_image,4);
else
    min_pixel_size = i_p.Results.min_independent_size;
    ad_zamir = find_ad_zamir(high_passed_image,threshed_image,min_pixel_size,'debug',i_p.Results.debug);
end

if(i_p.Results.status_messages), disp('Done finding adhesion regions'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove adhesions outside mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist('cell_mask','var'))
    for i = 1:max(ad_zamir(:))
        assert(any(any(ad_zamir == i)), 'Error: can''t find ad number %d', i);
        this_ad = zeros(size(ad_zamir));
        this_ad(ad_zamir == i) = 1;
        if (sum(sum(this_ad & cell_mask)) == 0)
            ad_zamir(ad_zamir == i) = 0;
        end
    end
    if(i_p.Results.status_messages), disp('Done removing adhesions outside the cell edge'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for too many adhesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_nums = unique(ad_zamir)';
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
if ((length(ad_nums) - 1) > i_p.Results.max_adhesion_count)
    system(['touch ', fullfile(output_dir, 'Found_too_many_adhesions')]);
    error(['Found more (',num2str(max(ad_zamir(:))),') adhesions than', ...
        ' max adhesion count (',num2str(i_p.Results.max_adhesion_count),').']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find and fill holes in single adhesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
props = regionprops(ad_zamir,'Area');
large_ad_nums = find([props.Area] >= 4);
for this_num = large_ad_nums
    %first make a binary image of the current adhesion and then run imfill
    %to fill any holes present
    this_ad = zeros(size(ad_zamir));
    this_ad(ad_zamir == this_num) = 1;
    filled_ad = imfill(this_ad,'holes');
    
    ad_zamir(logical(filled_ad)) = this_num;
    if (i_p.Results.debug && mod(this_num,50)==0)
        disp(['Done filling holes in ',num2str(this_num), '/', num2str(length(ad_nums))]);
    end
end
if(i_p.Results.status_messages), disp('Done filling adhesion holes'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renumber adhesions to be sequential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_nums = unique(ad_zamir);
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(ad_nums)
    ad_zamir(ad_zamir == ad_nums(i)) = i - 1;
end
if(i_p.Results.status_messages), disp('Done renumbering adhesion regions'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build adhesion perimeters image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_zamir_perim = zeros(size(ad_zamir));
for i = 1:max(ad_zamir(:))
    assert(any(any(ad_zamir == i)), 'Error: can''t find ad number %d', i);
    this_ad = zeros(size(ad_zamir));
    this_ad(ad_zamir == i) = 1;
    ad_zamir_perim(bwperim(this_ad)) = i;
end
if(i_p.Results.status_messages), disp('Done building adhesion perimeters'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imwrite(double(ad_zamir)/2^16,fullfile(output_dir, filenames.adhesions),'bitdepth',16);
imwrite(double(ad_zamir_perim)/2^16,fullfile(output_dir, filenames.adhesions_perim),'bitdepth',16);
imwrite(im2bw(ad_zamir,0),fullfile(output_dir, filenames.adhesions_binary));

highlighted_image = create_highlighted_image(focal_normed, im2bw(ad_zamir_perim,0), 'color_map',[1,1,0]);
if (exist('cell_mask','var'))
    highlighted_image = create_highlighted_image(highlighted_image, bwperim(cell_mask),'color_map',[1,0,0]);
end
imwrite(highlighted_image,fullfile(output_dir, 'highlights.png'));

if (i_p.Results.paper_figures)
    col_range = (find(sum(ad_zamir),1,'first')-5):(find(sum(ad_zamir),1,'last')+5);
    col_range = col_range(col_range > 0 & col_range < size(ad_zamir,2));
    row_range = (find(sum(ad_zamir,2),1,'first')-5):(find(sum(ad_zamir,2),1,'last')+5);
    row_range = row_range(row_range > 0 & row_range < size(ad_zamir,1));
    
    normed_hp_image = (high_passed_image - min(high_passed_image(:)))/range(high_passed_image(:));
    normed_hp_image = normed_hp_image(row_range,col_range);
    imwrite(normed_hp_image,fullfile(output_dir,'high_passed_image.png'),'bitdepth',16);
    
    imwrite(highlighted_image(row_range,col_range,1:3),fullfile(output_dir,'highlights_cropped.png'));
    imwrite(focal_normed(row_range,col_range),fullfile(output_dir,'focal_cropped.png'));
end

if (nargout > 0)
    varargout{1} = struct('adhesions',im2bw(ad_zamir,0),'ad_zamir',ad_zamir);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshed_image = find_threshed_image(high_passed_image, filter_thresh, proximity_filter)

if (length(filter_thresh) == 1)
    threshed_image = high_passed_image >= filter_thresh;
else
    high_threshed_image = high_passed_image >= filter_thresh(2);
    high_threshed_image = remove_edge_adhesions(high_threshed_image);
    high_threshed_image = imdilate(high_threshed_image,strel('disk',proximity_filter));
    
    low_threshed_image = high_passed_image >= filter_thresh(1);
    low_thresh_bwlabel = bwlabel(low_threshed_image,4);
    
    overlap_labels = unique(low_thresh_bwlabel.*high_threshed_image);
    if (overlap_labels(1) == 0)
        overlap_labels = overlap_labels(2:end);
    end
    
    threshed_image = ismember(low_thresh_bwlabel,overlap_labels);
end
1;