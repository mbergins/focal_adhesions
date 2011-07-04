function find_adhesion_properties(focal_file,adhesions_file,varargin)
% FIND_ADHESION_PROPERTIES    deteremines and outputs the quantitative
%                             properties associated with the adhesions
%                             located in prior steps
%
%   find_adhesion_properties(ff,af,OPTIONS) determines the quantitative
%   properites of the adhesions identified in the file 'af', using
%   information in the focal image file 'ff', the properties are written
%   to a set of csv files located a folder named 'raw_data' in the same
%   directory as the focal image file 'ff'
%
%   Options:
%
%       -cell_mask: file which contains the cell mask, defaults to not
%        present
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_ADHESION_PROPERTIES';

i_p.addRequired('focal_file',@(x)exist(x,'file') == 2);
i_p.addRequired('adhesions_file',@(x)exist(x,'file') == 2);

i_p.parse(focal_file, adhesions_file);

i_p.addParamValue('output_dir', fileparts(focal_file), @(x)exist(x,'dir')==7);
i_p.addOptional('cell_mask',0,@(x)exist(x,'file') == 2);
i_p.addOptional('protrusion_file',0,@(x)exist(x,'file') == 2);
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(focal_file, adhesions_file, varargin{:});

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask',i_p.UsingDefaults)))
    cell_mask = imread(i_p.Results.cell_mask);
end

%read in and normalize the input focal adhesion image
focal_image = double(imread(focal_file));

%read in the labeled adhesions
adhesions = imread(adhesions_file);

%check if protrusion_file is specified, read it in if specified
if (isempty(strmatch('protrusion_file',i_p.UsingDefaults)))
    load(i_p.Results.protrusion_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('cell_mask','var'))
    if (exist('protrusion','var'))
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'cell_mask',cell_mask,'protrusion_data',protrusion,'debug',i_p.Results.debug);
    else
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'cell_mask',cell_mask,'debug',i_p.Results.debug);
    end
else
    if (exist('protrusion_matrix','var'))
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'protrusion_data',protrusion,'debug',i_p.Results.debug);
    else
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'debug',i_p.Results.debug);
    end
end
if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adhesion_props = collect_adhesion_properties(orig_I,labeled_adhesions,varargin)
% COLLECT_ADHESION_PROPERTIES    using the identified adhesions, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                adhesions
%
%   ad_p = collect_adhesion_properties(ad_I,c_m,orig_I) collects the
%   properties of the adhesions identified in the binary image 'ad_I',
%   using the cell mask in 'c_m' and the original focal image data in
%   'orig_I', returning a structure 'ad_p' containing properties
%
%   Properties Collected:
%       -all of the properties collected by regioprops(...,'all')
%       -the distance of each adhesion's centroid from the nearest cell
%        edge
%       -the average and variance of the normalized fluorescence signal
%        within each adhesion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'COLLECT_ADHESION_PROPERTIES';

i_p.addRequired('orig_I',@isnumeric);
i_p.addRequired('labeled_adhesions',@(x)isnumeric(x));

i_p.addParamValue('cell_mask',0,@(x)isnumeric(x) || islogical(x));
i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addParamValue('protrusion_data',0,@(x)iscell(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(labeled_adhesions,orig_I,varargin{:});

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask',i_p.UsingDefaults)))
    cell_mask = i_p.Results.cell_mask;
end

if (isempty(strmatch('protrusion_data',i_p.UsingDefaults)))
    protrusion_data = i_p.Results.protrusion_data;
end

adhesion_props = regionprops(labeled_adhesions,'Area','Centroid', ... 
    'Eccentricity','MajorAxisLength','MinorAxisLength','Orientation', ... 
    'PixelIdxList','Solidity');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centroid_data = [adhesion_props.Centroid];
centroid_x = centroid_data(1:2:length(adhesion_props)*2);
centroid_y = centroid_data(2:2:length(adhesion_props)*2);

ad_centroid = [mean(centroid_x),mean(centroid_y)];
adhesion_props(1).Adhesion_centroid = ad_centroid;

for i=1:max(labeled_adhesions(:))
    adhesion_props(i).Average_adhesion_signal = mean(orig_I(labeled_adhesions == i));
    adhesion_props(i).Variance_adhesion_signal = var(orig_I(labeled_adhesions == i));
    adhesion_props(i).Max_adhesion_signal = max(orig_I(labeled_adhesions == i));
    adhesion_props(i).Min_adhesion_signal = min(orig_I(labeled_adhesions == i));
    
    this_ad = labeled_adhesions;
    this_ad(labeled_adhesions ~= i) = 0;
    this_ad = logical(this_ad);
    background_region = logical(imdilate(this_ad,strel('disk',i_p.Results.background_border_size,0)));
    background_region = and(background_region,not(labeled_adhesions));
    if (exist('cell_mask','var'))
        background_region = and(background_region,cell_mask);
    end
    adhesion_props(i).Background_adhesion_signal = mean(orig_I(background_region));
    adhesion_props(i).Background_area = sum(background_region(:));
    adhesion_props(i).Background_corrected_signal = adhesion_props(i).Average_adhesion_signal - adhesion_props(i).Background_adhesion_signal;
    
    shrunk_region = logical(imerode(this_ad,strel('disk',1,0)));
    if (sum(shrunk_region(:)) == 0), shrunk_region = this_ad; end
    adhesion_props(i).Shrunk_area = sum(shrunk_region(:));
    adhesion_props(i).Shrunk_adhesion_signal = mean(orig_I(shrunk_region));
    adhesion_props(i).Shrunk_corrected_signal = adhesion_props(i).Shrunk_adhesion_signal - adhesion_props(i).Background_adhesion_signal;
    
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(labeled_adhesions(:)))]); end
end

adhesion_mask = im2bw(labeled_adhesions,0);
adhesion_props(1).Adhesion_mean_intensity = sum(sum(orig_I(adhesion_mask)))/sum(sum(adhesion_mask));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Extracted If Protrusion Data Available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('protrusion_data','var'))
    all_data = [];
    for i=1:size(protrusion_data,2)
        all_data = [all_data, sqrt(protrusion_data{i}(:,3).^2 + protrusion_data{i}(:,4).^2)']; %#ok<AGROW>
    end
    median_velo = median(all_data); %#ok<NASGU>
    
    for i=1:size(protrusion_data,2)
        protrusion_matrix = protrusion_data{i};
        for j=1:max(labeled_adhesions(:))
            dists = sqrt((protrusion_matrix(:,1) - adhesion_props(j).Centroid(1)).^2 + (protrusion_matrix(:,2) - adhesion_props(j).Centroid(2)).^2);
            sorted_dists = sort(dists);
            best_line_nums = find(dists <= sorted_dists(5), 5,'first');
            
            edge_projection = [];
            edge_speed = [];
            for k=1:length(best_line_nums)
                this_line_num = best_line_nums(k);
                
                adhesion_to_edge = [protrusion_matrix(this_line_num,1) - adhesion_props(j).Centroid(1), protrusion_matrix(this_line_num,2) - adhesion_props(j).Centroid(2)];
                adhesion_to_edge = adhesion_to_edge / sqrt(adhesion_to_edge(1)^2 + adhesion_to_edge(2)^2);
                edge_vector = protrusion_matrix(this_line_num,3:4);
                if (sqrt(edge_vector(1)^2 + edge_vector(2)^2) > 10) 
                    edge_vector = (edge_vector / sqrt(edge_vector(1)^2 + edge_vector(2)^2))  * 10;
                end
                
                edge_projection(k) = sqrt(sum(edge_vector.^2))*(dot(edge_vector,adhesion_to_edge)/(sqrt(sum(edge_vector.^2)) * sqrt(sum(adhesion_to_edge.^2)))); %#ok<AGROW>
                edge_speed(k) = sqrt(sum(edge_vector.^2)); %#ok<AGROW>
            end
            adhesion_props(j).Edge_projection(i,1) = mean(edge_projection);
            adhesion_props(j).Edge_speed(i,1) = mean(edge_speed);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Extracted If Cell Mask Available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('cell_mask','var'))
    cell_centroid = regionprops(bwlabel(cell_mask),'centroid');
    cell_centroid = cell_centroid.Centroid;
    
    [border_row,border_col] = ind2sub(size(cell_mask),find(bwperim(cell_mask)));
    adhesion_props(1).Border_pix = [border_col,border_row];
    
    adhesion_props(1).Cell_size = sum(cell_mask(:));
    
    adhesion_props(1).Cell_mean_intensity = sum(sum(orig_I(cell_mask)))/adhesion_props(1).Cell_size;
    
    cell_not_ad_mask = cell_mask & not(adhesion_mask);
    adhesion_props(1).Cell_not_ad_mean_intensity = sum(sum(orig_I(cell_not_ad_mask)))/sum(sum(cell_not_ad_mask));
    
    not_cell_mask = not(cell_mask);
    adhesion_props(1).Outside_mean_intensity = sum(sum(orig_I(not_cell_mask)))/sum(sum(not_cell_mask));
    
    [dists, indices] = bwdist(~cell_mask);
    
    %Now we search for the pixels which are closest to an edge of the cell
    %mask that is also touching the edge of image. We want to find these
    %pixels because the true closest cell edge may be off the microscope
    %field of view. To be safe, we will set those
    %distance-to-nearest-cell-edge values to NaN.
    black_border_mask = cell_mask;
    black_border_mask(1,:) = 0; black_border_mask(end,:) = 0; 
    black_border_mask(:,1) = 0; black_border_mask(:,end) = 0;
    
    [bb_dists, bb_indexes] = bwdist(~black_border_mask);
    
    dists(bb_dists < dists) = NaN;
    for i=1:max(labeled_adhesions(:))
        centroid_pos = round(adhesion_props(i).Centroid);
        centroid_unrounded = adhesion_props(i).Centroid;
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_adhesion_properties - centroid not found');
            adhesion_props(i).Centroid_dist_from_edge = NaN;
        else
            adhesion_props(i).Centroid_dist_from_edge = dists(centroid_pos(2),centroid_pos(1));
            [cep_x,cep_y] = ind2sub(size(cell_mask), indices(centroid_pos(2), centroid_pos(1)));
            adhesion_props(i).Closest_edge_pixel = [cep_x,cep_y];
            
            adhesion_props(i).Centroid_dist_from_center = sqrt((cell_centroid(1) - centroid_unrounded(1))^2 + (cell_centroid(2) - centroid_unrounded(2))^2);
            adhesion_props(i).Angle_to_center = acos((centroid_unrounded(1) - cell_centroid(1))/adhesion_props(i).Centroid_dist_from_center);
            assert(adhesion_props(i).Angle_to_center >= 0 && adhesion_props(i).Angle_to_center <= pi, 'Error: angle to center out of range: %d',adhesion_props(i).Angle_to_center);
            if (centroid_unrounded(2) - cell_centroid(2) < 0)
                if (centroid_unrounded(1) - cell_centroid(1) < 0)
                    assert(adhesion_props(i).Angle_to_center >= pi/2 && adhesion_props(i).Angle_to_center <= pi)
                    adhesion_props(i).Angle_to_center = 2*pi - adhesion_props(i).Angle_to_center;
                elseif (centroid_unrounded(1) - cell_centroid(1) >= 0)
                    assert(adhesion_props(i).Angle_to_center >= 0 && adhesion_props(i).Angle_to_center <= pi/2)
                    adhesion_props(i).Angle_to_center = 2*pi - adhesion_props(i).Angle_to_center;
                end
            end
        end
        adhesion_props(i).CB_corrected_signal = adhesion_props(i).Average_adhesion_signal - adhesion_props(1).Cell_not_ad_mean_intensity;
    end
end

function write_adhesion_data(S,varargin)
% WRITE_STRUCT_DATA     write most the data stored in a given struct to a
%                       set of ascii formated files, using the field
%                       names as the file names
%
%   write_struct_data(S) writes most the field names in struct 'S' to ascii
%   formated files, suitable for use in other programs, the fieldnames are
%   used as the file names, all of the files are placed in subfolder
%   'raw_data' of the current working directory
%
%   write_struct_data(S,'out_dir',d) writes most the field names in struct
%   'S' to ascii formated files, suitable for use in other programs, the
%   fieldnames are used as the file names, all of the files are placed in
%   the provided directory 'd'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'WRITE_ADHESION_DATA';

i_p.addRequired('S',@isstruct);
i_p.addParamValue('out_dir','raw_data',@ischar);

i_p.parse(S,varargin{:});
out_dir = i_p.Results.out_dir;

if (not(exist(out_dir,'dir')))
    mkdir(out_dir);
end

to_exclude = {'ConvexHull','ConvexImage','Image','FilledImage', ...
    'PixelList', 'SubarrayIdx', 'Border_pix', 'Extrema'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names = fieldnames(S);

ad_props_cell = struct2cell(S);

print_strings = struct('PixelIdxList','%0.f');

for i = 1:size(field_names,1)
    
    if(strmatch(field_names(i),to_exclude))
        continue;
    end
    
    format_string = '%f';
    if(strmatch(field_names(i),fieldnames(print_strings)))
        format_string = print_strings.(field_names{i});
    end    
    
    file_out = fullfile(out_dir,[cell2mat(field_names(i)),'.csv']);
    
    data = ad_props_cell(i,:);
    output_CSV_from_cell(data, file_out, 'format', format_string);
end


function output_CSV_from_cell(data, out_file, varargin)
% output_CSV_from_cell    writes a provided cell data structure to a CSV
%                         file
%
%   write_struct_data(D,OF) writes the data in cell 'D' to file 'DF', using
%   fprintf format '%f'
%
%   Options:
%       -format: specify the string passed to fprintf for number conversion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'output_CSV_from_cell';

i_p.addRequired('data',@iscell);
i_p.addRequired('out_file', @(x) exist(fileparts(x),'dir') == 7 );

i_p.addParamValue('format','%f',@ischar);

i_p.parse(data,out_file,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_handle = fopen(i_p.Results.out_file,'wt');
for i = 1:max(size(data))
    for j = 1:max(size(data{i}))
        assert(any(size(data{i}) == 1))
        if (j < max(size(data{i})))
            fprintf(file_handle,[i_p.Results.format,','],data{i}(j));
        else
            fprintf(file_handle,[i_p.Results.format,'\n'],data{i}(j));
        end        
    end
end
fclose(file_handle);
