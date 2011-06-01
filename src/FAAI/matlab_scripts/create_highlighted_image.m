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

i_p.addParamValue('color_map',jet(double(max(high(:)))),@(x)(all(high(:) == 0) || (isnumeric(x) && (size(x,1) >= max(unique(high))))));
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
