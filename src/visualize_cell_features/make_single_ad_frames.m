function make_single_ad_frames(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('min_longevity',-Inf,@(x)isnumeric(x) && x > 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('..'))

filenames = add_filenames_to_struct(struct());

image_padding_min = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

individual_images_dir = fullfile(exp_dir,filenames.individual_results_dir);
image_folders = dir(individual_images_dir);
image_folders = image_folders(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect General Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracking_seq = csvread(fullfile(exp_dir,filenames.tracking)) + 1;
tracking_seq_size = size(tracking_seq);

% if (isempty(strmatch('adhesion_file', i_p.UsingDefaults)))
%     ad_to_include = csvread(i_p.Results.adhesion_file);
%     
%     temp_tracking_mat = zeros(tracking_seq_size);
%     temp_tracking_mat(ad_to_include(:,1),:) = tracking_seq(ad_to_include(:,1),:);
%     tracking_seq = temp_tracking_mat;
% end

longevities = sum(tracking_seq > 0,2);
tracking_seq(longevities < i_p.Results.min_longevity,:) = 0;

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

i_size = [];

for i_num = 1:length(image_folders)
    ad_label = imread(fullfile(individual_images_dir,image_folders(i_num).name,filenames.adhesions));
    
    if (i_num == 1)
        i_size = size(ad_label);
    end
    
    bounds = regionprops(ad_label,'BoundingBox');

    for i = 1:size(tracking_seq,1)
        tracking_row = tracking_seq(i,:);
        if (tracking_row(i_num) <= 0), continue; end
                
        ad_num = tracking_row(i_num);

        corners = [bounds(ad_num).BoundingBox(1), bounds(ad_num).BoundingBox(2)];
        corners = [corners, corners + bounds(ad_num).BoundingBox(3:4)]; %#ok<AGROW>

        if (corners(1) < bounding_matrix(i,1)), bounding_matrix(i,1) = corners(1); end
        if (corners(2) < bounding_matrix(i,2)), bounding_matrix(i,2) = corners(2); end
        if (corners(3) > bounding_matrix(i,3)), bounding_matrix(i,3) = corners(3); end
        if (corners(4) > bounding_matrix(i,4)), bounding_matrix(i,4) = corners(4); end
    end

%     if (mod(i_num,10) == 0 || j == max_image_num)
%         disp(['Bounding image: ',num2str(i_num)]);
%     end
end

bounding_matrix(:,1:2) = bounding_matrix(:,1:2) - image_padding_min;
bounding_matrix(:,1:2) = floor(bounding_matrix(:,1:2));
bounding_matrix(bounding_matrix(:,1) <= 0,1) = 1;
bounding_matrix(bounding_matrix(:,2) <= 0,2) = 1;

bounding_matrix(:,3:4) = bounding_matrix(:,3:4) + image_padding_min;
bounding_matrix(:,3:4) = ceil(bounding_matrix(:,3:4));
bounding_matrix(bounding_matrix(:,3) > i_size(2),3) = i_size(2);
bounding_matrix(bounding_matrix(:,4) > i_size(1),4) = i_size(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Single Ad Image Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_set_range = csvread(fullfile(individual_images_dir,image_folders(i_num).name,filenames.focal_image_min_max));

%all_images will hold all the highlighted frames produced for the movies,
%each row will hold the images from each lineages, where the columns will
%hold all the frames from each time point
all_images = cell(size(tracking_seq,1), 1);

for i_num = 1:length(image_folders)
    %Gather and scale the input adhesion image
    orig_i = double(imread(fullfile(individual_images_dir,image_folders(i_num).name,filenames.focal_image)));
    orig_i = (orig_i - image_set_range(1))/(image_set_range(2) - image_set_range(1));
    
    %Gather the ad label image
    ad_label_perim = imread(fullfile(individual_images_dir,image_folders(i_num).name,filenames.adhesions_perim));

    %Gather the cell edge image if available
    cell_mask_file = fullfile(individual_images_dir,image_folders(i_num).name,filenames.cell_mask);
    if (exist(cell_mask_file,'file'))
        cell_edge = bwperim(imread(cell_mask_file));
    end
    
    %cycle through each adhesion that we want to make a small multiple of
    for row_num = 1:size(tracking_seq,1)
        padded_num = sprintf(['%0',num2str(length(num2str(tracking_seq_size(1)))),'d'],row_num);
        
        tracking_row = tracking_seq(row_num,:);

        %now we do a check to see if there is an adhesion in the next
        %image or the image before, because we also want to render the
        %image immediately before birth and right after death
        surrounding_entries = [0, tracking_row(i_num), 0];
        try surrounding_entries(1) = tracking_row(i_num - 1); end %#ok<TRYNC>
        try surrounding_entries(3) = tracking_row(i_num + 1); end %#ok<TRYNC>
        
        if (any(surrounding_entries > 0))
            ad_num = tracking_row(i_num);
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
            
            all_images{row_num}{i_num} = highlighted_image;
        end
        
        %check if this adhesion's small multiple set is complete
        if (all(surrounding_entries <= 0) || i_num == length(image_folders))
            
            %if the current image set is empty, we haven't hit the adhesion
            %yet in the image set, so continue to the next adhesion
            if (size(all_images{row_num}, 2) == 0), continue; end
            
            if (exist('ad_to_include','var') && size(ad_to_include,2) == 2)
                offset_row_num = find(ad_to_include(:,1) == row_num);
                image_counts = ad_to_include(offset_row_num,2); %#ok<FNDSB>
                if (isempty(all_images{row_num}{1}))
                    image_counts = image_counts + 1;
                elseif (tracking_row(i_num) <= 0)
                    image_counts = image_counts + 1;
                end
                
                [~, ad_file_name] = fileparts(i_p.Results.adhesion_file);
                if (not(isempty(regexpi(ad_file_name,'disassembly'))))
                    offset_type = 'disassembly';
                elseif (not(isempty(regexpi(ad_file_name,'assembly'))))
                    offset_type = 'assembly';
                else
                    warning('FA:fileName','Expected to find either disassembly or assembly in adhesion_file parameter.')
                end
                
                output_file = fullfile(exp_dir,'visualizations', offset_type, [padded_num, '.png']);
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
            
            output_file = fullfile(exp_dir,'visualizations', 'overall', [padded_num, '.png']);
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
    if (mod(i_num,10) == 0)
        disp(['Highlight image: ',num2str(i_num)]);
    end
end

toc;