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

i_p = inputParser;
i_p.FunctionName = 'FIND_FOCAL_ADHESIONS';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.parse(I_file);

i_p.addParamValue('cell_mask',0,@(x)exist(x,'file') == 2);
i_p.addParamValue('min_size',0.56,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('pixel_size',0.215051,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('filter_thresh',0.1,@isnumeric);
i_p.addParamValue('scale_filter_thresh',0,@(x)islogical(x) || (isnumeric(x) && (x == 1 || x == 0)));
i_p.addParamValue('output_dir', fileparts(I_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('output_file', 'adhesions.png', @ischar);
i_p.addParamValue('output_file_perim', 'adhesions_perim.png', @ischar);
i_p.addParamValue('output_file_binary', 'adhesions_binary.png', @ischar);
i_p.addParamValue('no_ad_splitting', 0, @(x) islogical(x) || x == 1 || x == 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask', i_p.UsingDefaults)))
    cell_mask = imread(i_p.Results.cell_mask);
end

%read in and normalize the input focal adhesion image
focal_image  = imread(I_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

if (i_p.Results.scale_filter_thresh)
    filter_thresh = i_p.Results.filter_thresh * (max(focal_image(:)) - min(focal_image(:)));
else
    filter_thresh = i_p.Results.filter_thresh;
end

%Add the folder with all the scripts used in this master program
addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_filt = fspecial('disk',i_p.Results.filter_size);
blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
high_passed_image = focal_image - blurred_image;
threshed_image = logical(im2bw(high_passed_image,filter_thresh));

%identify and remove adhesions on the immediate edge of the image
threshed_image = remove_edge_adhesions(threshed_image);

min_pixel_size = floor((sqrt(i_p.Results.min_size)/i_p.Results.pixel_size)^2);
if (i_p.Results.no_ad_splitting)
    %if splitting is off, there is no need to use the fancy watershed based
    %segmentation methods, just identify the connected areas
    ad_zamir = bwlabel(threshed_image,4);
else
    ad_zamir = find_ad_zamir(high_passed_image,threshed_image,min_pixel_size,'debug',i_p.Results.debug);
end
disp('Done finding adhesion regions')

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
% Find and fill holes in single adhesions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
disp('Done filling adhesion holes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renumber adhesions to be sequential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_nums = unique(ad_zamir);
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(ad_nums)
    ad_zamir(ad_zamir == ad_nums(i)) = i - 1;
end

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

addpath(genpath('..'))

scaled_image = focal_image;
scaled_image = scaled_image - min(focal_image(:));
scaled_image = scaled_image .* (1/max(scaled_image(:)));

highlighted_image = create_highlighted_image(scaled_image, im2bw(ad_zamir_perim,0));
if (exist('cell_mask','var'))
    highlighted_image = create_highlighted_image(highlighted_image, bwperim(cell_mask),'color_map',[1,0,0]);
end
imwrite(highlighted_image,fullfile(i_p.Results.output_dir, 'highlights.png'));

if (nargout > 0)
    varargout{1} = struct('adhesions',im2bw(ad_zamir,0),'ad_zamir',ad_zamir);
end
