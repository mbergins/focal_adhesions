function make_ghost_frames(cfg_file,varargin)
%MAKE_GHOST_FRAMES    Builds single ghost images of the adhesions in a
%                     given experiment
%
%   make_ghost_frames(cfg_file,options) builds ghost images from raw
%   experimental data, where files are placed and the movie config options
%   are coded in cfg_file
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'MAKE_GHOST_FRAMES';

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

tracking_seq = load(tracking_seq_file) + 1;

ghost_frames_count = Inf;

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

max_live_adhesions = find_max_live_adhesions(tracking_seq);

lineage_cmap = jet(max_live_adhesions);
lineage_to_cmap = zeros(size(tracking_seq,1),1);

time_cmap = jet(size(tracking_seq,2));
birth_time_to_cmap = zeros(size(tracking_seq,1),1);

i_seen = 0;

for i = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],i);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));
    ad_label_perim = imread(fullfile(I_folder,padded_i_num,adhesions_perim_filename));
    ad_nums = tracking_seq(tracking_seq(:,i_seen) > 0,i_seen);
    temp_ad_label = zeros(size(ad_label));
    temp_ad_label_perim = zeros(size(ad_label_perim));
    for j = 1:length(ad_nums)
        this_num = ad_nums(j);
        assert(any(any(ad_label == this_num)), 'Error: can''t find ad num %d in image number %d.',this_num,padded_i_num)

        temp_ad_label(ad_label == this_num) = this_num;
        temp_ad_label_perim(ad_label_perim == this_num) = this_num;
    end
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
            catch
                assert(~(any(taken_dists == max(taken_dists),1,'first')), 'Error: could not find a possible color number in image number %d',padded_i_num);
            end
        end
    end

    %Make sure all the live adhesions have had a number assigned to their
    %lineage
    assert(all(lineage_to_cmap(tracking_seq(:,i_seen) > 0) > 0), 'Error in assigning unique color codes');
    assert(all(birth_time_to_cmap(tracking_seq(:,i_seen) > 0) > 0), 'Error in assigning birth time color codes');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Adhesion Ghost Images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Save the label matrices
    if (exist('labels','var'))
        frame_size_count = size(labels,2);
        if (frame_size_count >= ghost_frames_count), frame_size_count = ghost_frames_count - 1; end
        for j = frame_size_count:-1:1
            labels(j+1).ad_perim = labels(j).ad_perim; %#ok<AGROW>
            labels(j+1).ad_filled = labels(j).ad_filled; %#ok<AGROW>
        end
        labels(1).ad_perim = ad_label_perim;
        labels(1).ad_filled = ad_label;
    else
        labels(1).ad_perim = ad_label_perim;
        labels(1).ad_filled = ad_label;
    end

    %Draw the ghost images
    if (i_seen == size(tracking_seq,2))
        highlighted_ghost_unique = ones(i_size);
        highlighted_ghost_time = ones(i_size);

        highlighted_ghost_unique_filled = ones(i_size);
        highlighted_ghost_time_filled = ones(i_size);
        
        if(i_p.Results.debug), disp('Building Images'); end
        
        for m=size(labels,2):-1:1
            this_i_num = i_seen - m + 1;
            this_ad_perim = labels(m).ad_perim;
            this_ad_filled = labels(m).ad_filled;

            these_ad_nums = tracking_seq(tracking_seq(:,this_i_num) > 0,this_i_num);

            mix_percent = (size(labels,2) - m + 1)/size(labels,2);
            mix_percent = 1;

            %Unique colored adhesion lineage image drawing
            cmap_nums = lineage_to_cmap(tracking_seq(:,this_i_num) > 0);
            assert(length(these_ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in unique color ghost image creation %d',this_i_num);
            this_cmap = zeros(max(this_ad_perim(:)),3);
            this_cmap(these_ad_nums,:) = lineage_cmap(cmap_nums,:);
            highlighted_ghost_unique = create_highlighted_image(highlighted_ghost_unique,this_ad_perim,'color_map',this_cmap,'mix_percent',mix_percent);
            highlighted_ghost_unique_filled = create_highlighted_image(highlighted_ghost_unique_filled,this_ad_filled,'color_map',this_cmap,'mix_percent',mix_percent);

            %Birth time colored adhesion lineage image drawing
            cmap_nums = birth_time_to_cmap(tracking_seq(:,this_i_num) > 0);
            assert(length(these_ad_nums) == length(cmap_nums),'Error: the number of adhesions does not match the number of lineage numbers in birth time color ghost image creation %d',this_i_num);
            this_cmap = zeros(max(this_ad_perim(:)),3);
            this_cmap(these_ad_nums,:) = time_cmap(cmap_nums,:);
            highlighted_ghost_time = create_highlighted_image(highlighted_ghost_time,this_ad_perim,'color_map',this_cmap,'mix_percent',mix_percent);
            highlighted_ghost_time_filled = create_highlighted_image(highlighted_ghost_time_filled,this_ad_filled,'color_map',this_cmap,'mix_percent',mix_percent);
            
            if(i_p.Results.debug), disp(this_i_num); end
        end
        
        ghost_image_dir = fullfile(out_path,'ghost_images');
        if (not(exist(ghost_image_dir,'dir')))
            mkdir(fullfile(out_path,'ghost_images')); 
        end
        
        %Apply the bounding box
        highlighted_ghost_unique = highlighted_ghost_unique(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        highlighted_ghost_time = highlighted_ghost_time(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        highlighted_ghost_unique_filled = highlighted_ghost_unique_filled(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        highlighted_ghost_time_filled = highlighted_ghost_time_filled(b_box(2):b_box(4), b_box(1):b_box(3), 1:3);
        
        %Write the images out before making adding the scale bar
        imwrite(highlighted_ghost_unique,fullfile(ghost_image_dir,'ghost_uni_ns.png'));
        imwrite(highlighted_ghost_time,fullfile(ghost_image_dir,'ghost_time_ns.png'));
        imwrite(highlighted_ghost_unique_filled,fullfile(ghost_image_dir,'ghost_uni_filled_ns.png'));
        imwrite(highlighted_ghost_time_filled,fullfile(ghost_image_dir,'ghost_time_filled_ns.png'));
        
        if (exist('pixel_size','var'))
            highlighted_ghost_unique = draw_scale_bar(highlighted_ghost_unique,pixel_size);
            highlighted_ghost_time = draw_scale_bar(highlighted_ghost_time,pixel_size);
            highlighted_ghost_unique_filled = draw_scale_bar(highlighted_ghost_unique_filled,pixel_size);
            highlighted_ghost_time_filled = draw_scale_bar(highlighted_ghost_time_filled,pixel_size);
            imwrite(highlighted_ghost_unique,fullfile(ghost_image_dir,'ghost_uni.png'));
            imwrite(highlighted_ghost_time,fullfile(ghost_image_dir,'ghost_time.png'));
            imwrite(highlighted_ghost_unique_filled,fullfile(ghost_image_dir,'ghost_uni_filled.png'));
            imwrite(highlighted_ghost_time_filled,fullfile(ghost_image_dir,'ghost_time_filled.png'));            
        end
    end
    
    if(i_p.Results.debug), disp(i_seen); end
end

profile off;
if (i_p.Results.debug), profile viewer; end
