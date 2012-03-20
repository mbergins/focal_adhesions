function find_cell_mask(I_file)
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

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);

i_p.parse(I_file);

mask_image = double(imread(I_file));

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

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
threshed_mask = mask_image > intensity(imin(min_index));

%%Mask Cleanup
connected_areas = bwlabel(threshed_mask);%
region_sizes = regionprops(connected_areas, 'Area');

%filter out connected regions smaller than 2.5% of the total image area
threshed_mask = ismember(connected_areas, find([region_sizes.Area] > 0.025*length(mask_image(:))));

threshed_mask = imfill(threshed_mask,'holes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_file = fullfile(fileparts(I_file),filenames.cell_mask);

imwrite(threshed_mask, out_file)