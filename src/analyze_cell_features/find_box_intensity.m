function find_box_intensity(cfg_file,varargin)
%MAKE_SINGLE_AD_FRAMES    Builds single image montages that track single
%                         adhesions through their entire lifecycle,
%                         including frames immediately preceding and
%                         following the adhesion's lifetime, if available
%
%   make_single_ad_frames(cfg_file,options) builds single adhesion montages
%   from raw experimental data, most config options are set in cfg_file
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'MAKE_SINGLE_AD_FRAMES';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,varargin{:});
if (i_p.Results.debug == 1), profile off; profile on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(cfg_file);
while 1
    line = fgetl(fid);
    if ~ischar(line), break; end
    eval(line);
end

addpath(genpath(path_folders));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect General Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_image_num = find_max_image_num(I_folder);
folder_char_length = length(num2str(max_image_num));
i_size = size(imread(fullfile(I_folder,num2str(max_image_num),focal_image)));

%Add one to the tracking matrix to conpensate for matrix indexing beginning
%at zero in perl
tracking_seq = load(tracking_seq_file) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather Bounding Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The bounding matrix will hold the proper ranges to cut out only the
%relavent area for each adhesion's life cycle. The matrix is initialized
%with values that will always be replaced and the data is structured like
%this:
%   Column 1: top corner left (row)
%   Column 2: top corner left (column)
%   Column 3: bottom corner right (row)
%   Column 4: bottom corner right (column)
bounding_matrix = [Inf*ones(size(tracking_seq,1),1), Inf*ones(size(tracking_seq,1),1), ...
    -Inf*ones(size(tracking_seq,1),1), -Inf*ones(size(tracking_seq,1),1)];

%i_seen will keep track of the number of images that have actually been
%read into the program, we need to keep track of this due to skipped
%frames, which will show up as missing files in the following loop
i_seen = 0;

for j = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],j);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;

    ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));

    bounds = regionprops(ad_label,'BoundingBox');

    for i = 1:size(tracking_seq,1)
        tracking_row = tracking_seq(i,:);
        if (tracking_row(i_seen) <= 0), continue; end

        ad_num = tracking_row(i_seen);

        corners = [bounds(ad_num).BoundingBox(1), bounds(ad_num).BoundingBox(2)];
        corners = [corners, corners + bounds(ad_num).BoundingBox(3:4)]; %#ok<AGROW>

        if (corners(1) < bounding_matrix(i,1)), bounding_matrix(i,1) = corners(1); end
        if (corners(2) < bounding_matrix(i,2)), bounding_matrix(i,2) = corners(2); end
        if (corners(3) > bounding_matrix(i,3)), bounding_matrix(i,3) = corners(3); end
        if (corners(4) > bounding_matrix(i,4)), bounding_matrix(i,4) = corners(4); end
    end

    if (i_p.Results.debug && (mod(i_seen,10) == 0 || j == max_image_num))
        disp(['Bounding image: ',num2str(i_seen)]);
    end
end

bounding_matrix(:,1:2) = floor(bounding_matrix(:,1:2));
bounding_matrix(bounding_matrix(:,1) <= 0,1) = 1;
bounding_matrix(bounding_matrix(:,2) <= 0,2) = 1;

bounding_matrix(:,3:4) = ceil(bounding_matrix(:,3:4));
bounding_matrix(bounding_matrix(:,3) > i_size(2),3) = i_size(2);
bounding_matrix(bounding_matrix(:,4) > i_size(1),4) = i_size(1);

assert(all(all(isnumeric(bounding_matrix))),'Error: part of bounding matrix is not numeric')
assert(all(all(bounding_matrix >= 1)),'Error: part of bounding matrix is less than 1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Box Intensity Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%i_seen will keep track of the number of images that have actually been
%read into the program, we need to keep track of this due to skipped
%frames, which will show up as missing files in the following loop
i_seen = 0;

%intensity_values will hold the intensity boxed intensity values from each
%single adhesion
intensity_values = ones(size(tracking_seq)) * NaN;
for j = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],j);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;

    %Gather and scale the input adhesion image
    orig_i = imread(fullfile(I_folder,padded_i_num,focal_image));
    scale_factor = double(intmax(class(orig_i)));
    orig_i = double(orig_i)/scale_factor;

    %Each row of the track matrix holds the adhesion number sequence that
    %defines an adhesion's life cycle, cycle through this for the current
    %image number, extracting the box intensity at each point
    for i = 1:size(tracking_seq,1)
        
        %we skip these entries because they indicate dead or yet to be
        %born adhesions
        if (tracking_seq(i,i_seen) < 1), continue; end
        
        this_b_box = bounding_matrix(i,:);
        
        intensity_values(i,i_seen) = mean(mean(orig_i(this_b_box(2):this_b_box(4), this_b_box(1):this_b_box(3))));
    end
    if (i_p.Results.debug && (mod(i_seen,10) == 0 || j == max_image_num))
        disp(['Bounding image: ',num2str(i_seen)]);
    end
end

if (not(exist(lin_time_series_folder,'dir'))), mkdir(lin_time_series_folder); end
csvwrite(fullfile(lin_time_series_folder,'Box_intensity.csv'),intensity_values);

profile off;
if (i_p.Results.debug), profile viewer; end
