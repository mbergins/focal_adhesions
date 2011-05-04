function image_with_bar = draw_scale_bar(image_no_bar,pixel_size,varargin)
%DRAW_SCALE_BAR    Adds a scale bar to an image
%
%   draw_scale_bar(image,p_size,options) draws a scale bar on the provided
%   matlab variable, 'image', using the physical pixel size , 'p_size', in
%   microns to appropriately size the bar, various settings can be changed
%   with options:
%
%   Options:
%       -bar_size: default 10 microns, desired bar_size in microns
%       -position_code: default 2, lower right hand corner, desired scale
%        bar position in the image, codes follow:
%
%   Position Codes:
%       -1: Upper right hand corner
%       -2: Lower right hand corner
%       -3: Lower left hand corner
%       -4: Upper left hand corner
%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'DRAW_SCALE_BAR';

i_p.addRequired('image_no_bar',@(x)isnumeric(x));
i_p.addRequired('pixel_size',@(x)isnumeric(x) && x > 0);

i_p.addParamValue('bar_size',10,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('position_code',2,@(x)isnumeric(x) && x > 0 && x < 4);
i_p.addParamValue('bar_color',1,@(x)isnumeric(x) && x >= 0 && x <=1 );
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(image_no_bar,pixel_size,varargin{:});

%Pull out the command line variables
debug = i_p.Results.debug;
bar_size = i_p.Results.bar_size;
position_code = round(i_p.Results.position_code);

offset_standard = 0.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[image_height image_width image_depth] = size(image_no_bar);

bar_width = round(bar_size/pixel_size);
bar_height = round(0.01*image_height);
if (bar_height < 2)
    bar_height = 2;
end

if (isa(image_no_bar,'double'))
    bar = ones(bar_height,bar_width)*i_p.Results.bar_color;
else 
    bar = double(intmax(class(image_no_bar)))*ones(bar_height,bar_width)*i_p.Results.bar_color;
end
bar_size = size(bar);

%Determine Position Code
best_position_code = 2;

if (not(position_code))
    position_code = best_position_code;
end

switch position_code
    case 1
        row_offset = offset_standard;
        col_offset = 1-offset_standard;
    case 2
        row_offset = 1-offset_standard;
        col_offset = 1-offset_standard;
    case 3
        row_offset = 1-offset_standard;
        col_offset = offset_standard;
    case 4
        row_offset = offset_standard;
        col_offset = offset_standard;
end

%Draw the image with the bar
image_with_bar = image_no_bar;

for i=1:image_depth
    switch position_code
        case 1
            bar_rows = ceil(row_offset*image_height):ceil(row_offset*image_height)-1+bar_size(1);
            bar_cols = ceil(col_offset*image_width)-bar_size(2):ceil(col_offset*image_width)-1;
        case 2
            bar_rows = ceil(row_offset*image_height)-bar_size(1):ceil(row_offset*image_height)-1;
            bar_cols = ceil(col_offset*image_width)-bar_size(2):ceil(col_offset*image_width)-1;
        case 3
            bar_rows = ceil(row_offset*image_height)-bar_size(1):ceil(row_offset*image_height)-1;
            bar_cols = ceil(col_offset*image_width):ceil(col_offset*image_width)-1+bar_size(2);
        case 4
            bar_rows = ceil(row_offset*image_height):ceil(row_offset*image_height)-1+bar_size(1);
            bar_cols = ceil(col_offset*image_width):ceil(col_offset*image_width)-1+bar_size(2);
        otherwise
            error('ERROR: draw_scale_bar - unexpected position code, must be between 0-4');
    end
    image_with_bar(bar_rows,bar_cols,i) = bar;
end

end
