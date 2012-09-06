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
i_p.StructExpand = true;
i_p.FunctionName = 'FIND_ADHESION_PROPERTIES';

i_p.addRequired('focal_file',@(x)exist(x,'file') == 2);
i_p.addRequired('adhesions_file',@(x)exist(x,'file') == 2);

i_p.addParamValue('output_dir', fileparts(focal_file), @(x)exist(x,'dir')==7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(focal_file, adhesions_file, varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%read in the cell mask image if present
if (exist(fullfile(fileparts(focal_file),filenames.cell_mask),'file'))
    cell_mask = imread(fullfile(fileparts(focal_file),filenames.cell_mask));
end

%read in the kinase file if available
if (exist(fullfile(fileparts(focal_file),filenames.kinase),'file'))
    kinase = imread(fullfile(fileparts(focal_file),filenames.kinase));
end

%read in the input focal adhesion image data
focal_image = double(imread(focal_file));
adhesions = imread(adhesions_file);

%gather the background correction, if available
background_correction = 0;

cor_file = fullfile(fileparts(focal_file),filenames.background_intensity);
if (exist(cor_file,'file'))
    disp('Found background correction file, using to correct adhesion intensities');
    background_correction = csvread(cor_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('cell_mask','var'))
    adhesion_properties = collect_adhesion_properties(focal_image,adhesions,background_correction, ...
        'cell_mask',cell_mask,'debug',i_p.Results.debug);
else
    adhesion_properties = collect_adhesion_properties(focal_image,adhesions,background_correction, ...
        'debug',i_p.Results.debug);
end

if (exist('kinase','var'))
    adhesion_properties = collect_kinase_properties(kinase,adhesions,adhesion_properties);
end

if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adhesion_properties = collect_kinase_properties(kinase,adhesions,adhesion_properties)

adhesion_properties(1).Kinase_mean_intensity = mean(kinase(adhesions > 0));
kinase_prop = regionprops(adhesions,kinase,'MeanIntensity');
[adhesion_properties.Kinase_intensity] = kinase_prop.MeanIntensity;
adhesion_properties(1).FA_Kinase_intensity_corr = ...
    corr([adhesion_properties.Average_adhesion_signal]',[kinase_prop.MeanIntensity]');

for i=1:max(adhesions(:))
    this_ad = adhesions == i;
    
    surrounding = imdilate(this_ad, strel('disk',4)) & not(adhesions);
    adhesion_properties(i).Kinase_intensity_corrected = adhesion_properties(i).Kinase_intensity - ...
        mean(kinase(surrounding));
end

function adhesion_props = collect_adhesion_properties(orig_I,labeled_adhesions, ...
    background_correction,varargin)
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
i_p.addRequired('background_correction',@(x)isnumeric(x));

i_p.addParamValue('cell_mask',0,@(x)isnumeric(x) || islogical(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(labeled_adhesions,orig_I,background_correction,varargin{:});

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask',i_p.UsingDefaults)))
    cell_mask = i_p.Results.cell_mask;
end

adhesion_props = regionprops(labeled_adhesions,'Area','Centroid', ...
    'Eccentricity','MajorAxisLength','MinorAxisLength','Orientation', ...
    'PixelIdxList');

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
dist_to_centroid = sqrt((centroid_x - ad_centroid(1)).^2 ...
    + (centroid_y - ad_centroid(2)).^2);

x_dist_to_cent = centroid_x - ad_centroid(1);
y_dist_to_cent = ad_centroid(2) - centroid_y;

angle_to_ad_cent = atan2(y_dist_to_cent,x_dist_to_cent)*(180/pi);

convex_hull = bwconvhull(labeled_adhesions > 0);
convex_dists = bwdist(~convex_hull);

for i=1:max(labeled_adhesions(:))
    adhesion_props(i).Average_adhesion_signal = mean(orig_I(labeled_adhesions == i)) - background_correction;
    
    adhesion_props(i).Dist_to_FA_cent = dist_to_centroid(i);
    adhesion_props(i).Angle_to_FA_cent = angle_to_ad_cent(i);
    
    centroid_rounded = round(adhesion_props(i).Centroid);
    adhesion_props(i).CHull_dist = convex_dists(centroid_rounded(2),centroid_rounded(1));
    
    angle_to_FA_vector = [cosd(angle_to_ad_cent(i)),sind(angle_to_ad_cent(i))];
    FA_orientation_vector = [cosd(adhesion_props(i).Orientation),sind(adhesion_props(i).Orientation)];
    FA_orientation_vector_alt = [cosd(adhesion_props(i).Orientation + 180),
        sind(adhesion_props(i).Orientation + 180)];
    
    angle_between = acosd(dot(angle_to_FA_vector,FA_orientation_vector));
    angle_between_alt = acosd(dot(angle_to_FA_vector,FA_orientation_vector_alt));
    
    if (abs(angle_between) < abs(angle_between_alt))
        adhesion_props(i).Angle_diff_from_radial = angle_between;
    else
        adhesion_props(i).Angle_diff_from_radial = angle_between_alt;
    end
    
    if (i_p.Results.debug)
        temp = zeros(size(labeled_adhesions));
        temp(labeled_adhesions == i) = 1;
        temp(labeled_adhesions ~= i) = 2;
        temp(not(labeled_adhesions)) = 0;
        temp(round(ad_centroid(2))-3:round(ad_centroid(2))+3,round(ad_centroid(1))-3:round(ad_centroid(1))+3) = 3;
        
        if (mod(i,10) == 0)
            disp(['Finished Ad: ',num2str(i), '/', num2str(max(labeled_adhesions(:)))]);
        end
    end
end

adhesion_mask = im2bw(labeled_adhesions,0);
adhesion_props(1).Adhesion_mean_intensity = sum(sum(orig_I(adhesion_mask)))/sum(sum(adhesion_mask));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Extracted If Cell Mask Available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist('cell_mask','var'))
    cell_centroid = regionprops(cell_mask,'centroid');
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
    
    [bb_dists, ~] = bwdist(~black_border_mask);
    
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

print_strings = struct('PixelIdxList','%0.f','Angle_diff_from_radial','%0.2f',...
    'Orientation','%0.2f','Angle_to_FA_cent','%0.2f','MajorAxisLength','%0.2f',...
    'MinorAxisLength','%0.2f','Average_adhesion_signal','%0.2f', ...
    'Variance_adhesion_signal','%0.2f','Min_adhesion_signal','%0.2f', ...
    'Max_adhesion_signal','%0.2f','Kinase_intensity','%0.2f', ...
    'Kinase_intensity_corrected','%0.2f');

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
