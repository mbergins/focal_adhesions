function adhesion_matrix = make_ad_matrix(ad_size, intensity, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'make_ad_matrix';

i_p.addRequired('ad_size',@(x)isnumeric(x) && length(x) == 2);
i_p.addRequired('intensity',@isnumeric);

i_p.parse(ad_size,intensity,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%interp2 blows up on things smaller than 2, so short circuit at the
%beginnign for those types of adhesions
if (any(ad_size <= 2)) 
    adhesion_matrix = ones(ad_size)*intensity;
    return;
end

circle = zeros(max(i_p.Results.ad_size));

if (mod(size(circle,1), 2) == 0)
    points = [size(circle,1)/2, size(circle,1)/2+1];
    circle(points(1), points(1)) = 1;
    circle(points(1), points(2)) = 1;
    circle(points(2), points(1)) = 1;
    circle(points(2), points(2)) = 1;
else
    circle(size(circle,1)/2+0.5, size(circle,1)/2+0.5) = 1;
end

[X,Y] = meshgrid(1:size(circle,1));
[XI,YI] = meshgrid(1:(size(circle,1)-1)/(ad_size(2)-1):size(circle,1), ... 
                   1:(size(circle,1)-1)/(ad_size(1)-1):size(circle,1));

interped_circle = interp2(X,Y,bwdist(circle),XI,YI);

edge_point = [size(interped_circle,1), size(interped_circle,2)/2];
if (mod(edge_point(1),1) ~= 0)
    edge_point(1) = edge_point(1) + 0.5;
end

if (mod(edge_point(2),1) ~= 0)
    edge_point(2) = edge_point(2) + 0.5;
end

adhesion_matrix = intensity*(interped_circle <= interped_circle(edge_point(1),edge_point(2)));