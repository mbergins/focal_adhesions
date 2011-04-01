function bbox = find_time_series_bbox(folder,varargin)
% FIND_TIME_SERIES_BBOX    determines the minumum bounding box for a given
%                          image time series
%
%   [bbox] = find_time_series_bbox(f) searches through the folders in the
%   provided folder 'f' for folders containing files named 'adhesions.png',
%   these files are used to calculate the overall bounding box for an image
%   time series
%
%   [bbox] = find_time_series_bbox(f,'image_filename',im) searches through
%   the folders in the provided folder 'f' for folders containing files
%   named 'im', these files are used to calculate the overall bounding box
%   for an image time series


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('folder',@(x)exist(x,'dir')==7);
i_p.addParamValue('image_filename','cell_mask.png',@ischar);
i_p.addParamValue('backup_image_filename','adhesions.png',@ischar);

i_p.parse(folder,varargin{:});

image_filename = i_p.Results.image_filename;

folders_to_exclude = {'.','..'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_files = dir(folder);
num_files = size(all_files,1);

bbox = [Inf, Inf, -Inf, -Inf];

for i = 1:num_files
    full_folder_path = fullfile(folder,all_files(i).name);
    if (not(exist(full_folder_path,'dir')))
        continue;
    end
    if (strmatch(all_files(i).name,folders_to_exclude))
        continue;
    end
    
    best_image_path = fullfile(full_folder_path,i_p.Results.image_filename);
    
    if (not(exist(best_image_path,'file')))
        best_image_path = fullfile(full_folder_path,i_p.Results.backup_image_filename);
        
        if (not(exist(best_image_path,'file')))
            continue;
        end
    end
    
    this_bbox = find_binary_bounding_box(imread(best_image_path));
    
    if (bbox(1) > this_bbox(1)), bbox(1) = this_bbox(1); end
    if (bbox(2) > this_bbox(2)), bbox(2) = this_bbox(2); end
    if (bbox(3) < this_bbox(3)), bbox(3) = this_bbox(3); end
    if (bbox(4) < this_bbox(4)), bbox(4) = this_bbox(4); end
    
end
