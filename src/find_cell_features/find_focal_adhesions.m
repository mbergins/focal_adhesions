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
%       -filter_thresh: threshold used to identify focal adhesions in the
%        average filtered image, defaults to 0.1
%       -output_dir: folder used to hold all the results, defaults to the
%        same folder as the image file, 'I'
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'FIND_FOCAL_ADHESIONS';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.parse(I_file);

i_p.addParamValue('cell_mask',0,@(x)exist(x,'file') == 2);

%Adhesion filtering parameters
i_p.addParamValue('min_adhesion_size',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('filter_thresh',0.1,@isnumeric);
i_p.addParamValue('min_independent_size',14,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('no_ad_splitting', 0, @(x) islogical(x) || x == 1 || x == 0);
i_p.addParamValue('max_adhesion_count', Inf, @(x) isnumeric(x));

%output parameters
i_p.addParamValue('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('output_file', 'adhesions.png', @ischar);
i_p.addParamValue('output_file_perim', 'adhesions_perim.png', @ischar);
i_p.addParamValue('output_file_binary', 'adhesions_binary.png', @ischar);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');
addpath(genpath('..'));

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask', i_p.UsingDefaults)))
    cell_mask = imread(i_p.Results.cell_mask);
end

filenames = add_filenames_to_struct(struct());

filter_thresh = csvread(fullfile(fileparts(I_file),filenames.focal_image_threshold));

%read in and normalize the input focal adhesion image
focal_image  = double(imread(I_file));
image_set_min_max = csvread(fullfile(fileparts(I_file),filenames.focal_image_min_max));
focal_normed = (focal_image - image_set_min_max(1))/(range(image_set_min_max));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply filter to find adhesion regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_filt = fspecial('disk',i_p.Results.filter_size);
blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
high_passed_image = focal_image - blurred_image;
threshed_image = high_passed_image > filter_thresh;

%identify and remove adhesions on the immediate edge of the image
threshed_image = remove_edge_adhesions(threshed_image);

%filter out small adhesions if requested
if (i_p.Results.min_adhesion_size > 1)
    labeled_thresh = bwlabel(threshed_image,4);
    
    props = regionprops(labeled_thresh,'Area'); %#ok<MRPBW>
    labeled_thresh = ismember(labeled_thresh, find([props.Area] >= i_p.Results.min_adhesion_size));
    
    threshed_image = labeled_thresh > 0;
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
if (i_p.Results.no_ad_splitting || not(isempty(strmatch('min_independent_size', i_p.UsingDefaults))))
    %if splitting is off, there is no need to use the fancy watershed based
    %segmentation methods, just identify the connected areas
    ad_zamir = bwlabel(threshed_image,4);
else
    min_pixel_size = i_p.Results.min_independent_size;
    ad_zamir = find_ad_zamir(high_passed_image,threshed_image,min_pixel_size,'debug',i_p.Results.debug);
end

disp('Done finding adhesion regions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renumber adhesions to be sequential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_nums = unique(ad_zamir);
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(ad_nums)
    ad_zamir(ad_zamir == ad_nums(i)) = i - 1;
end
disp('Done renumbering adhesion regions')

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
    disp('Done removing adhesions outside the cell edge')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for too many adhesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (max(ad_zamir(:)) > i_p.Results.max_adhesion_count)
    system(['touch ', fullfile(i_p.Results.output_dir, 'Found_too_many_adhesions')]);
    error(['Found more (',num2str(max(ad_zamir(:))),') adhesions than', ...
        ' max adhesion count (',num2str(i_p.Results.max_adhesion_count),').']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find and fill holes in single adhesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
ad_nums = unique(ad_zamir);
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(ad_nums)
    this_num = ad_nums(i);
    
    %first make a binary image of the current adhesion and then run imfill
    %to fill any holes present
    this_ad = zeros(size(ad_zamir));
    this_ad(ad_zamir == this_num) = 1;
    filled_ad = imfill(this_ad);
    
    ad_zamir(logical(filled_ad)) = this_num;
    if (i_p.Results.debug && mod(i,50)==0)
        disp(['Done filling holes in ',num2str(i), '/', num2str(length(ad_nums))]);
    end
end
toc;
disp('Done filling adhesion holes')

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
disp('Done building adhesion perimeters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(double(ad_zamir)/2^16,fullfile(i_p.Results.output_dir, i_p.Results.output_file),'bitdepth',16);
imwrite(double(ad_zamir_perim)/2^16,fullfile(i_p.Results.output_dir, i_p.Results.output_file_perim),'bitdepth',16);
imwrite(im2bw(ad_zamir,0),fullfile(i_p.Results.output_dir, i_p.Results.output_file_binary));

highlighted_image = create_highlighted_image(focal_normed, im2bw(ad_zamir_perim,0));
if (exist('cell_mask','var'))
    highlighted_image = create_highlighted_image(highlighted_image, bwperim(cell_mask),'color_map',[1,0,0]);
end
imwrite(highlighted_image,fullfile(i_p.Results.output_dir, 'highlights.png'));

if (nargout > 0)
    varargout{1} = struct('adhesions',im2bw(ad_zamir,0),'ad_zamir',ad_zamir);
end
toc;
