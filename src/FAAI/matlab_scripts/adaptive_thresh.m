function threshold = adaptive_thresh(I,varargin)
% ADPATIVE_THRESH    finds a threshold which separates a given image into
%                    background and signal, assuming a bimodal histrogram
%                    is present, uses the algoritm from p.599-600 in
%                    Digital Image Processing by Gonzalez and Woods
%
%   adaptive_thresh(I) using the pixel values in image 'I', find and return
%   a threshold which should separate the image into background and signal 
%
%   adaptive_thresh(I,'up_mean_weight',W) using the pixel values in image
%   'I', find and return a threshold which should separate the image into
%   background and signal using a fractional weight 'W' on the upper mean


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'ADAPTIVE_THRESH';

i_p.addRequired('I', @isnumeric);

i_p.addParamValue('upper_mean_weight',0.5,@(x) isnumeric(x) && x <= 1 && x >= 0);

i_p.parse(I,varargin{:});

upper_mean_weight = i_p.Results.upper_mean_weight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = mean(I(:));

threshold_change_thresh = 0.0001;

upper_mean = mean(I(find(I > threshold)));
lower_mean = mean(I(find(I <= threshold)));

new_thresh = upper_mean_weight*upper_mean + (1-upper_mean_weight)*lower_mean;

while (abs(threshold - new_thresh) > threshold_change_thresh)
    threshold = new_thresh;
    
    upper_mean = mean(I(find(I > threshold)));
    lower_mean = mean(I(find(I <= threshold)));

    new_thresh = upper_mean_weight*upper_mean + (1-upper_mean_weight)*lower_mean;
end

end