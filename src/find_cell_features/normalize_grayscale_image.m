function [I] = normalize_grayscale_image(I,varargin)
% NORMALIZE_GRAYSCALE_IMAGE    using a given image, the same image with the
%                              minimum and maximum pixel value adjusted to
%                              0 and 1 respectively is returned in double
%                              format
%
%   normalize_grayscale_image(I) normalize the image 'I' to have min and
%   max values of 0 and 1
%
%   normalize_grayscale_image(I,'min_max',M) normalize the image 'I'
%   assuming that the two values in the array 'M' are the min and max
%   value, with the first value 'M(1)' is the min and the second 'M(2)' is
%   the max

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'NORMALIZE_GRAYSCALE_IMAGE';

i_p.addRequired('I',@isnumeric);

i_p.parse(I);

i_p.addParamValue('min_max',[min(I(:)),max(I(:))],@(x) size(x,2)==2 && isnumeric(x(1)) && isnumeric(x(2)));

i_p.parse(I,varargin{:});

min_max = i_p.Results.min_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_size = size(I);
I = double(I);
min_max = double(min_max);

if (size(image_size,2) > 2)
    if (image_size(3) > 1)
        I = I(:,:,1);
    end
end

I = I - min_max(1);
I = I/(min_max(2)-min_max(1));

end