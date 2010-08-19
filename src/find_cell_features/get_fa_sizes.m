function get_fa_sizes(I_file,varargin)
% GET_FA_SIZES    locates the focal adhesions in a given image, optionally
%                 returns the segmented image or writes the segmented image
%                 to a file
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
i_p.FunctionName = 'GET_FA_SIZES';

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.parse(I_file);

i_p.addParamValue('pixel_size',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('filter_thresh',0.1,@isnumeric);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_file,varargin{:});

I_filt = fspecial('disk',i_p.Results.filter_size);

profile off;
if (i_p.Results.debug)
    profile -timer real;
    profile on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_num = length(imfinfo(I_file));

min_max = [Inf, -Inf];
for i=1:image_num
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in and normalize the input focal adhesion image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    focal_image  = imread(I_file,i);
    scale_factor = double(intmax(class(focal_image)));
    focal_image  = double(focal_image)/scale_factor;
    
    this_min = min(focal_image(:));
    this_max = max(focal_image(:));
    if (min_max(1) > this_min)
        min_max(1) = this_min;
    end
    
    if (min_max(2) < this_max)
        min_max(2) = this_max;
    end    
end

%scale the filter thresh depending on the min and max values found in the
%images
filter_thresh = i_p.Results.filter_thresh *(min_max(2) - min_max(1));

fa_highlight_movie = avifile(fullfile(fileparts(I_file),'mask_highlight.avi'));
all_areas = [];
for i=1:image_num
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in and normalize the input focal adhesion image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    focal_image  = imread(I_file,i);
    scale_factor = double(intmax(class(focal_image)));
    focal_image  = double(focal_image)/scale_factor;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the cell mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    cell_mask = find_cell_mask(focal_image);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter and threshold the FA image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
    high_passed_image = focal_image - blurred_image;
    threshed_image = logical(im2bw(high_passed_image,filter_thresh));
    
    %identify and remove adhesions on the immediate edge of the image
    threshed_image = remove_edge_adhesions(threshed_image);
    
    ad = bwlabel(threshed_image,4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove adhesions outside mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:max(ad(:))
        assert(any(any(ad == j)), 'Error: can''t find ad number %d', j);
        this_ad = zeros(size(ad));
        this_ad(ad == j) = 1;
        if (sum(sum(this_ad & cell_mask)) == 0)
            ad(ad == j) = 0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find and fill holes in single adhesions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_nums = unique(ad);
    assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
    for j = 2:length(ad_nums)
        this_num = ad_nums(j);
        
        %first make a binary image of the current adhesion and then run imfill
        %to fill any holes present
        this_ad = zeros(size(ad));
        this_ad(ad == this_num) = 1;
        filled_ad = imfill(this_ad);
        
        ad(logical(filled_ad)) = this_num;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Renumber adhesions to be sequential
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_nums = unique(ad);
    assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
    for j = 2:length(ad_nums)
        ad(ad == ad_nums(j)) = j - 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build adhesion perimeters image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_perim = zeros(size(ad));
    for j = 1:max(ad(:))
        assert(any(any(ad == j)), 'Error: can''t find ad number %d', j);
        this_ad = zeros(size(ad));
        this_ad(ad == j) = 1;
        ad_perim(bwperim(this_ad)) = j;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add frame to the diagnostic movie
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    scaled_image = (focal_image - min_max(1)) / (min_max(2) - min_max(1));
    
    highlighted_image = create_highlighted_image(scaled_image, im2bw(ad_perim,0));
    highlighted_image = create_highlighted_image(highlighted_image, bwperim(cell_mask),'color_map',[1,0,0]);
    fa_highlight_movie = addframe(fa_highlight_movie, highlighted_image);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collect and Output FA Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_props = regionprops(ad,'Area');
    
    all_areas = [all_areas, [temp_props.Area]*i_p.Results.pixel_size^2]; %#ok<AGROW>
    
    if (mod(i,10) == 0)
        disp(['Done with Image ', num2str(i), '/', num2str(image_num)]);
    end
end

dlmwrite(fullfile(fileparts(I_file),'areas.csv'),all_areas);

fa_highlight_movie = close(fa_highlight_movie); %#ok<NASGU>

if (i_p.Results.debug)
    profile viewer;
end

toc
end

function cleaned_binary = remove_edge_adhesions(threshed_image)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'REMOVE_EDGE_ADHESIONS';

i_p.addRequired('threshed_image',@(x)isnumeric(x) || islogical(x));

i_p.parse(threshed_image);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edge_border = ones(size(threshed_image));
edge_border = bwperim(edge_border);

[row_indexes,col_indexes] = ind2sub(size(threshed_image), find(edge_border));
edge_adhesions = bwselect(threshed_image,col_indexes,row_indexes,4);

cleaned_binary = threshed_image & not(edge_adhesions);

end

function [varargout] = find_cell_mask(mask_image)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image
%
%   create_cell_mask_image(I,OF) finds the cell mask using the image in
%   file 'I' and writes the binary cell mask to the output file 'OF'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_CELL_MASK';

i_p.addRequired('mask_image',@(x)isnumeric(x));

i_p.parse(mask_image);

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Threshold identification
sorted_mask_pixels = sort(mask_image(:));
% sorted_mask_pixels(1:0.05*round(length(sorted_mask_pixels))) = 0;

[heights, intensity] = hist(sorted_mask_pixels,1000);

smoothed_heights = smooth(heights,0.05,'loess');
[zmax,imax,zmin,imin]= extrema(smoothed_heights);

%keep in mind that the zmax is sorted by value, so the highest peak is
%first and the corresponding index is also first in imax, the same pattern
%hold for zmin and imin

sorted_max_indexes = sort(imax);
first_max_index = find(sorted_max_indexes == imax(1));

%locate the index between the first two maximums
min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));

threshed_mask = im2bw(mask_image, intensity(imin(min_index)));

%%Mask Cleanup
connected_areas = bwlabel(threshed_mask);%
region_sizes = regionprops(connected_areas, 'Area'); %#ok<MRPBW>

%filter out connected regions smaller than 10 pixels
threshed_mask = ismember(connected_areas, find([region_sizes.Area] == max([region_sizes.Area])));

threshed_mask = imfill(threshed_mask,'holes');

if (nargout >= 1)
    varargout{1} = threshed_mask;
end

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

i_p.addParamValue('color_map',[0,1,0],@(x)(all(high(:) == 0) || (isnumeric(x) && (size(x,1) >= max(unique(high))))));
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

%check if the highlight image is blank, if so create the highlight image,
%without highlights and return
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

function [xmax,imax,xmin,imin] = extrema(x)
%EXTREMA   Gets the global extrema points from a time series.
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA(X) returns the global minima and maxima 
%   points of the vector X ignoring NaN's, where
%    XMAX - maxima points in descending order
%    IMAX - indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - indexes of the XMIN
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema).
%
%   Example:
%      x = 2*pi*linspace(-1,1);
%      y = cos(x) - 0.5 + 0.5*rand(size(x)); y(40:45) = 1.85; y(50:53)=NaN;
%      [ymax,imax,ymin,imin] = extrema(y);
%      plot(x,y,x(imax),ymax,'g.',x(imin),ymin,'r.')
%
%   See also EXTREMA2, MAX, MIN

%   Written by
%   Lic. on Physics Carlos Adri�n Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico, 2004
%
%   nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish. 
% 2006-11-17 : Accept NaN's.
% 2007-04-09 : Change name to MAXIMA, and definition added.


xmax = [];
imax = [];
xmin = [];
imin = [];

% Vector input?
Nt = numel(x);
if Nt ~= length(x)
 error('Entry must be a vector.')
end

% NaN's:
inan = find(isnan(x));
indx = 1:Nt;
if ~isempty(inan)
 indx(inan) = [];
 x(inan) = [];
 Nt = length(x);
end

% Difference between subsequent elements:
dx = diff(x);

% Is an horizontal line?
if ~any(dx)
 return
end

% Flat peaks? Put the middle element:
a = find(dx~=0);              % Indexes where x changes
lm = find(diff(a)~=1) + 1;    % Indexes where a do not changes
d = a(lm) - a(lm-1);          % Number of elements in the flat peak
a(lm) = a(lm) - floor(d/2);   % Save middle elements
a(end+1) = Nt;

% Peaks?
xa  = x(a);             % Serie without flat peaks
b = (diff(xa) > 0);     % 1  =>  positive slopes (minima begin)  
                        % 0  =>  negative slopes (maxima begin)
xb  = diff(b);          % -1 =>  maxima indexes (but one) 
                        % +1 =>  minima indexes (but one)
imax = find(xb == -1) + 1; % maxima indexes
imin = find(xb == +1) + 1; % minima indexes
imax = a(imax);
imin = a(imin);

nmaxi = length(imax);
nmini = length(imin);                

% Maximum or minumim on a flat peak at the ends?
if (nmaxi==0) && (nmini==0)
 if x(1) > x(Nt)
  xmax = x(1);
  imax = indx(1);
  xmin = x(Nt);
  imin = indx(Nt);
 elseif x(1) < x(Nt)
  xmax = x(Nt);
  imax = indx(Nt);
  xmin = x(1);
  imin = indx(1);
 end
 return
end

% Maximum or minumim at the ends?
if (nmaxi==0) 
 imax(1:2) = [1 Nt];
elseif (nmini==0)
 imin(1:2) = [1 Nt];
else
 if imax(1) < imin(1)
  imin(2:nmini+1) = imin;
  imin(1) = 1;
 else
  imax(2:nmaxi+1) = imax;
  imax(1) = 1;
 end
 if imax(end) > imin(end)
  imin(end+1) = Nt;
 else
  imax(end+1) = Nt;
 end
end
xmax = x(imax);
xmin = x(imin);

% NaN's:
if ~isempty(inan)
 imax = indx(imax);
 imin = indx(imin);
end

% Same size as x:
imax = reshape(imax,size(xmax));
imin = reshape(imin,size(xmin));

% Descending order:
[temp,inmax] = sort(-xmax); clear temp
xmax = xmax(inmax);
imax = imax(inmax);
[xmin,inmin] = sort(xmin);
imin = imin(inmin);


% Carlos Adri�n Vargas Aguilera. nubeobscura@hotmail.com
end