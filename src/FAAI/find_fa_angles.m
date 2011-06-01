function find_fa_angles(exp_dir,varargin)
tstart=tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('stdev_thresh',2,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('min_axial_ratio',3,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('min_adhesion_size',3,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('color_blind',0,@(x)isnumeric(x) && x == 1 || x == 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Threshold Determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh_start = tic;
find_exp_thresholds(exp_dir,'stdev_thresh',i_p.Results.stdev_thresh);
disp('Done with finding threshold')
toc(thresh_start);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FA Finding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fa_finding_start = tic;
base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

for i_num = 1:size(image_dirs,1)
    puncta_file = fullfile(base_dir,image_dirs(i_num).name,filenames.focal_image);
    find_focal_adhesions(puncta_file,'status_messages',0,'min_adhesion_size', ...
        i_p.Results.min_adhesion_size);
    
    if (mod(i_num,1) == 0)
        fprintf('Done find FAs: %d/%d\n',i_num, size(image_dirs,1))
        time_left = round(toc(fa_finding_start)*((size(image_dirs,1) - i_num)/(60*i_num)));
        fprintf('About %d minutes left\n',time_left)
    end
end
toc(fa_finding_start);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Property Extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
props_start = tic;
all_props = struct();
for i_num = 1:size(image_dirs,1)
    image_set = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
    props = regionprops(image_set.adhesions,'MajorAxisLength','MinorAxisLength','Orientation','Area');
    
    propnames = fieldnames(props);
    for i = 1:length(propnames)
        this_prop = propnames{i};
        if (isempty(strmatch(this_prop,fieldnames(all_props))))
            all_props.(this_prop) = [];
        end
        all_props.(this_prop) = [all_props.(this_prop), props.(this_prop)];
        dlmwrite(fullfile(base_dir,image_dirs(i_num).name,[this_prop,'.csv']),[props.(this_prop)]);
    end
    
    ratio = [props.MajorAxisLength]./[props.MinorAxisLength];
    low_ratio = find(ratio < i_p.Results.min_axial_ratio);
    high_ratio = setdiff(1:max(image_set.adhesions(:)), low_ratio);
    
    low_ratio_perim = ismember(image_set.adhesions_perim,low_ratio);
    high_ratio_perim = ismember(image_set.adhesions_perim,high_ratio);
    
    highlight = image_set.focal_norm;
    if (i_p.Results.color_blind)
        highlight = create_highlighted_image(highlight,high_ratio_perim,'color_map',[1,1,0],'mix_percent',0.5);
        highlight = create_highlighted_image(highlight,low_ratio_perim,'color_map',[0,0,1],'mix_percent',0.5);        
    else
        highlight = create_highlighted_image(highlight,high_ratio_perim,'color_map',[0,1,0],'mix_percent',0.5);
        highlight = create_highlighted_image(highlight,low_ratio_perim,'color_map',[1,0,0],'mix_percent',0.5);
    end
    
    imwrite(highlight,fullfile(base_dir,image_dirs(i_num).name,'ratio_highlight.png'));
    
    if (mod(i_num,10) == 0)
        disp(sprintf('Done getting properites: %d/%d',i_num, size(image_dirs,1)))
    end
end

all_data_output_dir = fullfile(exp_dir,'adhesion_props','lin_time_series');
if (not(exist(all_data_output_dir,'dir'))) 
    mkdir(all_data_output_dir);
end

propnames = fieldnames(all_props);
for i = 1:length(propnames)
    this_prop = propnames{i};
    dlmwrite(fullfile(all_data_output_dir,[this_prop,'.csv']),[all_props.(this_prop)]);
end
toc(props_start);

toc(tstart);
1;