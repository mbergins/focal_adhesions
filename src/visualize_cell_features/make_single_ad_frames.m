function make_single_ad_frames(cfg_file,varargin)
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
tic;
i_p = inputParser;
i_p.FunctionName = 'MAKE_SINGLE_AD_FRAMES';

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('start_row',1,@(x)x >= 1);
i_p.addParamValue('end_row',1,@(x)x >= 1);
i_p.addParamValue('adhesion_file',@(x)exist(x,'file') == 2);

i_p.parse(cfg_file,varargin{:});

addpath(genpath('..'))
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
tracking_seq_size = size(tracking_seq);

%Filter the tracking matrix according to 'start_row' and 'end_row' rules
%provided via command parameters.
if (   isempty(strmatch('start_row', i_p.UsingDefaults)) ...
    && isempty(strmatch('end_row', i_p.UsingDefaults)))

    assert(i_p.Results.start_row <= size(tracking_seq,1));
    assert(i_p.Results.end_row <= size(tracking_seq,1));
    tracking_beginning = zeros(i_p.Results.start_row - 1,size(tracking_seq,2));
    tracking_middle = tracking_seq(i_p.Results.start_row:i_p.Results.end_row,:);
    tracking_end = zeros(size(tracking_seq,1) - i_p.Results.end_row, size(tracking_seq,2));
    
    tracking_seq = [tracking_beginning; tracking_middle; tracking_end];
    
elseif (isempty(strmatch('start_row', i_p.UsingDefaults)))

    assert(i_p.Results.start_row <= size(tracking_seq,1));
    tracking_beginning = zeros(i_p.Results.start_row - 1,size(tracking_seq,2));
    tracking_end = tracking_seq(i_p.Results.start_row:end,:);
    
    tracking_seq = [tracking_beginning; tracking_end];
    
elseif (isempty(strmatch('end_row', i_p.UsingDefaults)))
    
    assert(i_p.Results.end_row <= size(tracking_seq,1));
    tracking_beginning = tracking_seq(1:i_p.Results.end_row,:);
    tracking_end = zeros(size(tracking_seq,1) - i_p.Results.end_row, size(tracking_seq,2));
    
    tracking_seq = [tracking_beginning; tracking_end];
    
end
assert(all(size(tracking_seq) == tracking_seq_size));

if (isempty(strmatch('adhesion_file', i_p.UsingDefaults)))
    ad_to_include = csvread(i_p.Results.adhesion_file);
    
    temp_tracking_mat = zeros(tracking_seq_size);
    temp_tracking_mat(ad_to_include(:,1),:) = tracking_seq(ad_to_include(:,1),:);
    tracking_seq = temp_tracking_mat;
end

rows_to_examine = zeros(size(tracking_seq,1),1);
for i=1:size(tracking_seq)
    rows_to_examine(i) = any(tracking_seq(i,:) > 0);
end
rows_to_examine = find(rows_to_examine);

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

    if (mod(i_seen,10) == 0 || j == max_image_num)
        disp(['Bounding image: ',num2str(i_seen)]);
    end
end

bounding_matrix(:,1:2) = bounding_matrix(:,1:2) - single_image_padding_min;
bounding_matrix(:,1:2) = floor(bounding_matrix(:,1:2));
bounding_matrix(bounding_matrix(:,1) <= 0,1) = 1;
bounding_matrix(bounding_matrix(:,2) <= 0,2) = 1;

bounding_matrix(:,3:4) = bounding_matrix(:,3:4) + single_image_padding_min;
bounding_matrix(:,3:4) = ceil(bounding_matrix(:,3:4));
bounding_matrix(bounding_matrix(:,3) > i_size(2),3) = i_size(2);
bounding_matrix(bounding_matrix(:,4) > i_size(1),4) = i_size(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Single Ad Image Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_set_min_max = csvread(fullfile(I_folder,'01',filenames.focal_image_min_max));

%i_seen will keep track of the number of images that have actually been
%read into the program, we need to keep track of this due to skipped
%frames, which will show up as missing files in the following loop
i_seen = 0;

%all_images will hold all the highlighted frames produced for the movies,
%each row will hold the images from each lineages, where the columns will
%hold all the frames from each time point
all_images = cell(size(tracking_seq,1), 1);

for j = 1:max_image_num
    padded_i_num = sprintf(['%0',num2str(folder_char_length),'d'],j);

    if (not(exist(fullfile(I_folder,padded_i_num,focal_image),'file'))), continue; end

    i_seen = i_seen + 1;

    %Gather and scale the input adhesion image
    orig_i = double(imread(fullfile(I_folder,padded_i_num,focal_image)));
    orig_i = (orig_i - image_set_min_max(1))/range(image_set_min_max);
    
    %Gather the ad label image
    ad_label_perim = imread(fullfile(I_folder,padded_i_num,adhesions_perim_filename));

    %Gather the cell edge image if available
    if (exist(fullfile(I_folder,padded_i_num,edge_filename),'file'))
        cell_edge = bwperim(imread(fullfile(I_folder,padded_i_num,edge_filename)));
    end
    
    for i = 1:size(rows_to_examine,1)
        row_num = rows_to_examine(i);
        
        padded_num = sprintf(['%0',num2str(length(num2str(tracking_seq_size(1)))),'d'],row_num);
        
        tracking_row = tracking_seq(row_num,:);

        %now we do a check to see if there is an adhesion in the next
        %image or the image before, because we also want to render the
        %image immediately before birth and right after death
        surrounding_entries = [0, tracking_row(i_seen), 0];
        try surrounding_entries(1) = tracking_row(i_seen - 1); end %#ok<TRYNC>
        try surrounding_entries(3) = tracking_row(i_seen + 1); end %#ok<TRYNC>
        
        if (any(surrounding_entries > 0))
            ad_num = tracking_row(i_seen);
            if (ad_num <= 0); ad_num = -Inf; end
            bounded_ad_label_perim = ad_label_perim(bounding_matrix(row_num,2):bounding_matrix(row_num,4), ...
                bounding_matrix(row_num,1):bounding_matrix(row_num,3));
            
            this_ad = zeros(size(bounded_ad_label_perim));
            this_ad(bounded_ad_label_perim == ad_num) = 1;
            
            not_this_ad = xor(im2bw(bounded_ad_label_perim,0),this_ad);
            assert(sum(sum(not_this_ad)) + sum(sum(this_ad)) == sum(sum(im2bw(bounded_ad_label_perim,0))))
            
            highlighted_image = orig_i(bounding_matrix(row_num,2):bounding_matrix(row_num,4), bounding_matrix(row_num,1):bounding_matrix(row_num,3));
            highlighted_image = create_highlighted_image(highlighted_image,this_ad,'color_map',[0,1,0]);
            highlighted_image = create_highlighted_image(highlighted_image,not_this_ad,'color_map',[0,0,1]);
            if (exist('cell_edge','var'))
                bounded_edge = cell_edge(bounding_matrix(row_num,2):bounding_matrix(row_num,4), bounding_matrix(row_num,1):bounding_matrix(row_num,3));
                highlighted_image = create_highlighted_image(highlighted_image,bounded_edge,'color_map',[1,0,0]);
            end
            
            all_images{row_num}{i_seen} = highlighted_image;
        end
        
        if (all(surrounding_entries <= 0) || j == max_image_num)
            
            %exit out of this loop through the tracking matrix if all the
            %current image set is empty, we have yet to hit the adhesions
            if (size(all_images{row_num}, 2) == 0), continue; end
            
            if (exist('ad_to_include','var') && size(ad_to_include,2) == 2)
                offset_row_num = find(ad_to_include(:,1) == row_num);
                image_counts = ad_to_include(offset_row_num,2); %#ok<FNDSB>
                if (isempty(all_images{row_num}{1}))
                    image_counts = image_counts + 1;
                elseif (tracking_row(i_seen) <= 0)
                    image_counts = image_counts + 1;
                end
                
                [folder, ad_file_name] = fileparts(i_p.Results.adhesion_file);
                if (not(isempty(regexpi(ad_file_name,'disassembly'))))
                    offset_type = 'disassembly';
                elseif (not(isempty(regexpi(ad_file_name,'assembly'))))
                    offset_type = 'assembly';
                else
                    warning('FA:fileName','Expected to find either disassembly or assembly in adhesion_file parameter.')
                end
                
                
                output_file = fullfile(out_path_single, offset_type, [padded_num, '.png']);
                %Draw the scale bar
                if (exist('pixel_size','var'))
                    write_montage_image_set(all_images{row_num},output_file, ...
                        'phase',offset_type,'num_images',image_counts, ...
                        'pixel_size',pixel_size,'bar_size',5)
                else
                    write_montage_image_set(all_images{row_num},output_file, ...
                        'phase',offset_type,'num_images',image_counts) 
                end 
            end
            
            output_file = fullfile(out_path_single, 'overall', [padded_num, '.png']);
            if (exist('pixel_size','var'))
                write_montage_image_set(all_images{row_num},output_file, ...
                    'pixel_size',pixel_size,'bar_size',5);
            else
                write_montage_image_set(all_images{row_num},output_file);
            end
            
            all_images{row_num} = cell(0);
            continue;
        end
        
    end
    if (mod(j,10) == 0)
        disp(['Highlight image: ',num2str(i_seen)]);
    end
end

profile off;
if (i_p.Results.debug), profile viewer; end