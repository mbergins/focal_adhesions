function bbox = find_bounding_box(I)
% FIND_BOUNDING_BOX    finds the minimum box which completely encloses
%                      the provided image
%
%   [min_x,min_y,max_x,max_y] = find_binary_bounding_box(I) finds the min
%   and max, x and y coordinates that completely encloses the pixels in
%   image 'I', whose value are greater than zero


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('I',@(x)isnumeric(x) || islogical(x));

i_p.parse(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bbox = [Inf, Inf, -Inf, -Inf];

for i = 1:size(I,3)
    if (isnumeric(I))
        this_binary_layer = im2bw(I(:,:,i),0);
    else
        this_binary_layer = I;
    end
    
    this_bounding_set = find_binary_bounding_box(this_binary_layer);
    if (this_bounding_set(1) < bbox(1)), bbox(1) = this_bounding_set(1); end
    if (this_bounding_set(2) < bbox(2)), bbox(2) = this_bounding_set(2); end
    
    if (this_bounding_set(3) > bbox(3)), bbox(3) = this_bounding_set(3); end
    if (this_bounding_set(4) > bbox(4)), bbox(4) = this_bounding_set(4); end
end