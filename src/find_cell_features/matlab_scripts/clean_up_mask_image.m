function [mask_image, edge_binary_image] = clean_up_mask_image(input_edge_binary)
%CLEAN_UP_EDGE_IMAGE   isolates the largest enclosed area in a provided
%                       binary image and fills in any holes in the enclosed
%                       area and returns the edge of the filled in area
%
%   EBI = clean_up_mask_image(EB) finds the largest enclosed area in image
%   'EB', isolates that enclosed area, fills in the holes in that area and
%   returns the edge of that filled in area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'CLEAN_UP_MASK_IMAGE';

i_p.addRequired('input_edge_binary',@(x) isnumeric(x) || islogical(x));

i_p.parse(input_edge_binary);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connected_areas = bwlabel(input_edge_binary);
region_sizes = regionprops(connected_areas, 'Area');

%filter out connected regions smaller than 10 pixels
edge_binary_image = ismember(connected_areas, find([region_sizes.Area] > 10));

% edge_binary_image = imdilate(edge_binary_image,strel('diamond',1));

mask_image = imfill(edge_binary_image,'holes');

end
