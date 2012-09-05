function make_single_ad_folders(cfg_file,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,varargin{:});

if (i_p.Results.debug == 1), profile off; profile on; end

addpath('matlab_scripts');

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

max_live_adhesions = find_max_live_adhesions(tracking_seq);

lineage_cmap = jet(max_live_adhesions);
lineage_to_cmap = zeros(size(tracking_seq,1),1);

i_seen = 0;

for i = 1:max_image_num
    if(i_p.Results.debug), 
%         if (i_seen > 10)
%             continue;
%         end
    end
    if (i_seen + 1 > size(tracking_seq,2))
        continue;
    end
    
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],i);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;
    padded_i_seen = sprintf(['%0',num2str(folder_char_length),'d'],i_seen);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_label = imread(fullfile(I_folder,padded_i_num,adhesions_filename));
    ad_label_perim = imread(fullfile(I_folder,padded_i_num,adhesions_perim_filename));
    
    %DO NOT add the bounding box around this image, it messes up potrace,
    %crazily
%     ad_label = ad_label(b_box(2):b_box(4), b_box(1):b_box(3));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Build the matrices translating number to colormap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:size(tracking_seq,1)
        %if the adhesion idenfied by the current lineage is not alive, skip
        %this lineage
        if (tracking_seq(j,i_seen) <= 0), continue; end

        %Unique lineage colors
        if (lineage_to_cmap(j) == 0)
            used_c_nums = sort(lineage_to_cmap(tracking_seq(:,i_seen) > 0));
            used_c_nums = used_c_nums(used_c_nums ~= 0);

            taken_nums = zeros(1,max_live_adhesions);
            taken_nums(used_c_nums) = 1;
            taken_dists = bwdist(taken_nums);

            try
                lineage_to_cmap(j) = find(taken_dists == max(taken_dists),1,'first');
                
                scale_amount = i_seen/max_image_num;
                
                scaled_cmap(j,:) = scale_amount*lineage_cmap(lineage_to_cmap(j),:);
            catch
                assert(isempty(find(taken_dists == max(taken_dists),1,'first')), 'Error: could not find a possible color number in image number %d',padded_i_num);
            end
        end
    end

    
    if (not(exist(fullfile(out_path,'../single_ad_files'),'dir'))), mkdir(fullfile(out_path,'../single_ad_files')); end
    for j = 1:size(tracking_seq,1)
        %if the adhesion idenfied by the current lineage is not alive, skip
        %this lineage
        if (tracking_seq(j,i_seen) <= 0), continue; end
        
        this_ad = zeros(size(ad_label));
        
        this_ad(ad_label == tracking_seq(j,i_seen)) = 1;
        
        if (not(exist(fullfile(out_path,'../single_ad_files',num2str(j)),'dir'))) 
            mkdir(fullfile(out_path,'../single_ad_files',num2str(j))); 
        end
        
        imwrite(this_ad,fullfile(out_path,'../single_ad_files',num2str(j),[padded_i_seen, '.png']));
    end


    if(i_p.Results.debug), disp(i_seen); end
end

csvwrite(fullfile(out_path,'../single_ad_files','c_map.csv'), lineage_cmap(lineage_to_cmap,:));
csvwrite(fullfile(out_path,'../single_ad_files','scaled_c_map.csv'), scaled_cmap);

profile off;
if (i_p.Results.debug), profile viewer; end
