function bbox = find_binary_bounding_box(I)
% FIND_BINARY_BOUNDING_BOX    finds the minimum box which completely encloses
%                             the provided binary image
%
%   [min_row,min_col,max_row,max_col] = find_binary_bounding_box(I) finds the min
%   and max, x and y coordinates that completely enclose the binary 
%   image, 'I'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('I',@(x)isnumeric(x) || islogical(x));

i_p.parse(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bbox(1) = find(sum(I) > 0,1,'first');
bbox(2) = find(sum(I,2) > 0,1,'first');
bbox(3) = find(sum(I)>0, 1, 'last');
bbox(4) = find(sum(I,2)>0, 1, 'last');