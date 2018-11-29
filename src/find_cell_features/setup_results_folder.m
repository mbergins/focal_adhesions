function setup_results_folder(I_source,out_folder,out_name,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('I_source',@(x)exist(x,'dir')==7);
i_p.addRequired('out_folder',@(x)ischar(x));
i_p.addRequired('out_name',@(x)ischar(x));

i_p.addParamValue('debug',0,@(x)(isnumeric(x) && (x == 0 || x == 1) || islogical(x)));

i_p.parse(I_source,out_folder,out_name,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_names = dir(i_p.Results.I_source); i_names = i_names(3:end);

%Verify that multiple multi-frame images are present in the folder
if (length(i_names) > 1)
    for i_num = 1:length(i_names)
        i_full_file = fullfile(i_p.Results.I_source,i_names(i_num).name);
        
        %deal with the extra .DS_store files Mac OS X likes to leave around
        if (strcmp(i_names(i_num).name, '.DS_store'))
            next;
        end
        
        i_count = length(imfinfo(i_full_file));
        if (i_count > 1) 
            error(['Detected multiple files of which ' , i_full_file, ' contains multiple frames.']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copy Files into the Results Directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (not(exist(i_p.Results.out_folder,'dir'))), mkdir(i_p.Results.out_folder); end

i_total = 0;
for i_num = 1:length(i_names)
    i_full_file = fullfile(i_p.Results.I_source,i_names(i_num).name);
    i_count = length(imfinfo(i_full_file));
    for sub_i_num = 1:i_count
        i_total = i_total+1;
        
        if (i_count == 1) 
            input_image = uint16(imread(i_full_file));
        else
            input_image = uint16(imread(i_full_file,sub_i_num));
        end
        
        image_out_folder = fullfile(i_p.Results.out_folder,sprintf('%05d',i_total));
        if (not(exist(image_out_folder,'dir'))), mkdir(image_out_folder); end
        
        %check for images with more than 2 dimensions, in that case, cut it
        %down to just the first dimension
        if (length(size(input_image)) > 2)
            marked_channel = find_label_channel(input_image);
            input_image = input_image(:,:,marked_channel);
        end
        
        imwrite(input_image,fullfile(image_out_folder,i_p.Results.out_name),...
            'Bitdepth',16);
        
        temp = imread(fullfile(image_out_folder,i_p.Results.out_name));
        assert(all(all(input_image == temp)));
    end
end

toc;
end

function FA_channel = find_label_channel(input_image)

channel_sums = [sum(sum(input_image(:,:,:)))];
FA_channel = find(channel_sums == max(channel_sums));

end