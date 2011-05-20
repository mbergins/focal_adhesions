function get_fa_sizes(I_folder,varargin)
% GET_FA_SIZES    locates the focal adhesions in a set of images in a given
%                 folder, finds the size of each adhesion and the number of
%                 adhesions found per image
%                         
%
%   get_fa_sizes(dir,OPTIONS) locate the focal adhesions in the images in
%   directory, 'dir'
%
%   Options:
%       
%       -pixel_size: the size of the side of the pixels in the image in
%        microns, if specified output size values will be in square
%        microns, if not specified output size values are in pixels
%       -filter_size: size of the averaging filter to use, defaults to 23
%       -filter_thresh: threshold used to identify focal adhesions in the
%        high pass filtered image, defaults to 0.1
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'GET_FA_SIZES';

i_p.addRequired('I_folder',@(x)exist(x,'dir') == 7);

i_p.parse(I_folder);

i_p.addParamValue('pixel_size',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('filter_thresh',0.1,@isnumeric);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_folder,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_files = dir(I_folder);
assert(strcmp(image_files(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_files(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
image_files = image_files(3:end);

min_max = [Inf, -Inf];
for i=1:length(image_files)
    
    this_file = fullfile(I_folder,image_files(i).name);
    if (exist(this_file,'dir'))
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in and normalize the input focal adhesion image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    focal_image  = double(imread(this_file));
    
    this_min = min(focal_image(:));
    this_max = max(focal_image(:));
    if (min_max(1) > this_min)
        min_max(1) = this_min;
    end
    
    if (min_max(2) < this_max)
        min_max(2) = this_max;
    end
end

I_filt = fspecial('disk',i_p.Results.filter_size);

all_areas = [];
ad_counts = [];

for i=1:length(image_files)
    
    this_file = fullfile(I_folder,image_files(i).name);
    if (exist(this_file,'dir'))
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read in raw data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    focal_image = double(imread(this_file));
    focal_image = (focal_image - min_max(1))/range(min_max);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter and threshold the FA image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    blurred_image = imfilter(focal_image,I_filt,'same',mean(focal_image(:)));
    high_passed_image = focal_image - blurred_image;
    threshed_image = logical(im2bw(high_passed_image,i_p.Results.filter_thresh));
    
    %identify and remove adhesions on the immediate edge of the image
    threshed_image = remove_edge_adhesions(threshed_image);
    
    ad = bwlabel(threshed_image,4);
        
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
    % Output Image with Adhesions Marked
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c_map = jet(max(ad_perim(:)));
    c_map = c_map(randperm(max(ad_perim(:))),:);
    highlighted_image = create_highlighted_image(focal_image, ad_perim, 'color_map',c_map);
    
    [pathstr, filename] = fileparts(this_file);
    if (not(exist(fullfile(fileparts(this_file),'highlights'),'dir')))
        mkdir(fullfile(fileparts(this_file),'highlights'))
    end
    imwrite([cat(3,focal_image,focal_image,focal_image),highlighted_image], ...
        fullfile(fileparts(this_file),'highlights',[filename,'_highlights.png']));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collect and Output FA Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_props = regionprops(ad,'Area');
    
    all_areas = [all_areas, [temp_props.Area]*i_p.Results.pixel_size^2]; %#ok<AGROW>
    ad_counts = [ad_counts, max(ad(:))]; %#ok<AGROW>
    
    disp(['Done with Image: ',this_file]);    
end

if (not(exist(fullfile(I_folder,'data'),'dir')))
    mkdir(fullfile(I_folder,'data'))
end
dlmwrite(fullfile(I_folder,'data','areas.csv'),all_areas);
dlmwrite(fullfile(I_folder,'data','ad_counts.csv'),ad_counts);

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