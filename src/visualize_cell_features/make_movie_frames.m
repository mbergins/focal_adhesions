function make_movie_frames(cfg_file,varargin)
%MAKE_MOVIE_FRAMES    Builds movie frames with the adhesions highlighted in
%                     using various conventions
%
%   make_movie_frames(cfg_file,options) builds individual movie frames from
%   raw experimental data, where files are placed and the movie config
%   options are coded in cfg_file
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.FunctionName = 'MAKE_MOVIE_FRAMES';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('no_b_box',0,@(x) islogical(x) || x == 0 || x == 1);
i_p.addParamValue('no_scale_bar',0,@(x) islogical(x) || x == 0 || x == 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,varargin{:});

addpath(genpath('..'));

filenames = add_filenames_to_struct(struct());

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

tracking_seq = load(tracking_seq_file) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find edges of image data in adhesion images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist(bounding_box_file,'file'))
    b_box = load(bounding_box_file);
else
    b_box = find_time_series_bbox(I_folder);
    csvwrite(bounding_box_file,b_box);
end

b_box(1:2) = b_box(1:2) - image_padding_min;
b_box(3:4) = b_box(3:4) + image_padding_min;
if (b_box(1) <= 0), b_box(1) = 1; end
if (b_box(2) <= 0), b_box(2) = 1; end
if (b_box(3) > i_size(2)), b_box(3) = i_size(2); end
if (b_box(4) > i_size(1)), b_box(4) = i_size(1); end

if (i_p.Results.no_b_box)
    b_box(1) = 1;
    b_box(2) = 1;
    b_box(3) = i_size(2);
    b_box(4) = i_size(1);
end

edge_cmap = jet(size(tracking_seq,2));
%define the edge image here because the edge image will be added to each
%image loop, so the image should be global
edge_image_ad = ones(i_size(1),i_size(2),3);

max_live_adhesions = find_max_live_adhesions(tracking_seq);

lineage_cmap = jet(max_live_adhesions);
lineage_to_cmap = zeros(size(tracking_seq,1),1);

time_cmap = jet(size(tracking_seq,2));
birth_time_to_cmap = zeros(size(tracking_seq,1),1);

image_set_range = csvread(fullfile(I_folder,'01',filenames.focal_image_min_max));
i_seen = 0;

for i = 1:max_image_num
    if (i_seen + 1 > size(tracking_seq,2))
        continue;
    end
    
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],i);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;
    padded_i_seen = sprintf(['%0',num2str(folder_char_length),'d'],i_seen);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather and scale the input adhesion image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orig_i = double(imread(fullfile(I_folder,padded_i_num,focal_image)));
    orig_i = (orig_i - image_set_range(1))/(image_set_range(2) - image_set_range(1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));
    ad_label_perim = imread(fullfile(I_folder,padded_i_num,adhesions_perim_filename));
    ad_nums = tracking_seq(tracking_seq(:,i_seen) > 0,i_seen);
    temp_ad_label = zeros(size(ad_label));
    temp_ad_label_perim = zeros(size(ad_label_perim));
    temp_ad_label(ismember(ad_label,ad_nums)) = ad_label(ismember(ad_label,ad_nums));
    temp_ad_label_perim(ismember(ad_label,ad_nums)) = ad_label_perim(ismember(ad_label,ad_nums));

    ad_label = temp_ad_label;
    ad_label_perim = temp_ad_label_perim;
    
    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename))); %#ok<NASGU>
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Build the matrices translating number to colormap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:size(tracking_seq,1)
        %if the adhesion idenfied by the current lineage is not alive, skip
        %this lineage
        if (tracking_seq(j,i_seen) <= 0), continue; end

        %Add the current adhesion to the time dependent color map if no
        %number is currently defined
        if (birth_time_to_cmap(j) == 0), birth_time_to_cmap(j) = i_seen; end

        %Unique lineage colors
        if (lineage_to_cmap(j) == 0)
            used_c_nums = sort(lineage_to_cmap(tracking_seq(:,i_seen) > 0));
            used_c_nums = used_c_nums(used_c_nums ~= 0);

            taken_nums = zeros(1,max_live_adhesions);
            taken_nums(used_c_nums) = 1;
            taken_dists = bwdist(taken_nums);

            try
                lineage_to_cmap(j) = find(taken_dists == max(taken_dists),1,'first');
            catch map_error %#ok<NASGU>
                assert(~any(taken_dists == max(taken_dists)), 'Error: could not find a possible color number in image number %d',padded_i_num);
            end
        end
    end

    %Make sure all the live adhesions have had a number assigned to their
    %lineage
    assert(all(lineage_to_cmap(tracking_seq(:,i_seen) > 0) > 0), 'Error in assigning unique color codes');
    assert(all(birth_time_to_cmap(tracking_seq(:,i_seen) > 0) > 0), 'Error in assigning birth time color codes');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Build the unique lineage highlighted image
    cmap_nums = lineage_to_cmap(tracking_seq(:,i_seen) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in unique lineage numbers image %d',padded_i_num);
    this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums,:) = lineage_cmap(cmap_nums,:);
    highlighted_all = create_highlighted_image(orig_i,ad_label_perim,'color_map',this_cmap);

    %Build the birth time highlighted image
    cmap_nums = birth_time_to_cmap(tracking_seq(:,i_seen) > 0);
    assert(length(ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the color map indexes in birth time image %d',padded_i_num);
    this_cmap = zeros(max(ad_label_perim(:)),3);
    this_cmap(ad_nums,:) = time_cmap(cmap_nums,:);
    highlighted_time = create_highlighted_image(orig_i,ad_label_perim,'color_map',this_cmap);

    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
        edge_image_ad = create_highlighted_image(edge_image_ad,cell_edge,'color_map',edge_cmap(i_seen,:));
        highlighted_all = create_highlighted_image(highlighted_all,cell_edge,'color_map',edge_cmap(i_seen,:));
    end
    edge_image_ad = create_highlighted_image(edge_image_ad,im2bw(ad_label_perim,0),'color_map',edge_cmap(i_seen,:));
    
    %Bound the images according to the bounding box found towards the top
    %of the program
    orig_i = orig_i(b_box(2):b_box(4), b_box(1):b_box(3));
    highlighted_all = highlighted_all(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
    highlighted_time = highlighted_time(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
    edge_image_ad_bounded = edge_image_ad(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);

    spacer = 0.5*ones(size(edge_image_ad_bounded,1),1,3);
    
    frame = cell(1,3);
    frame{1} = [edge_image_ad_bounded,spacer,highlighted_all];
    frame{2} = [cat(3,orig_i,orig_i,orig_i),spacer,highlighted_all];
    frame{3} = [edge_image_ad_bounded,spacer,highlighted_time];
    
    %Add scale bars if the pixel size is available
    if (exist('pixel_size','var') && not(i_p.Results.no_scale_bar))
        for j = 1:size(frame,2)
            frame{j} = draw_scale_bar(frame{j},pixel_size);
        end
    end
    
    %Output the original unaltered image if requested
    if (output_original_image)
        if (not(exist(fullfile(out_path,'orig_i'),'dir'))), mkdir(fullfile(out_path,'orig_i')); end
        imwrite(orig_i,fullfile(out_path,'orig_i',[padded_i_seen,'.png']));
    end
    
    %Output all the other images
    if (exist('out_path','var'))
        for j = 1:length(out_prefix) %#ok<USENS>
            output_filename = fullfile(out_path,out_prefix{1,j},[padded_i_seen,'.png']);
            fullpath = fileparts(output_filename);
            if (not(exist(fullpath,'dir')))
                mkdir(fullpath);
            end
            imwrite(frame{j},output_filename);
        end
    end

    disp(['Done with ', num2str(i_seen), '/',num2str(max_image_num)]);
end

profile off;
if (i_p.Results.debug), profile viewer; end
