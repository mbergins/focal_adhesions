function mask_thresh = find_cell_mask(I_file,varargin)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image
%
%   create_cell_mask_image(I) finds the cell mask using the image in
%   file 'I' and writes the binary cell mask to the same folder as 'I'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.StructExpand = true;

i_p.addRequired('I_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('mask_threshold',0,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('median_filter',0,@(x)x==1 || x==0);

i_p.addParamValue('debug',0,@(x)x==1 || x==0);

i_p.parse(I_file,varargin{:});

mask_image = double(imread(I_file));

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (i_p.Results.median_filter)
    mask_image = medfilt2(mask_image,[3,3]);
end

mask_pixels = mask_image(:);
mask_pixels = mask_pixels(mask_pixels >= quantile(mask_pixels,0.05));

if (not(any(strcmp(i_p.UsingDefaults,'mask_threshold'))))
    mask_thresh = i_p.Results.mask_threshold;
else
    %%Threshold identification
    [heights, intensity] = hist(mask_pixels(:),length(unique(mask_pixels(:)))/4);
    
    smoothed_heights = smooth(heights,0.05,'loess');
    [zmax,imax,~,imin]= extrema(smoothed_heights);
    
    if (length(imax) == 1)
        mask_thresh = NaN;
    else
        %keep in mind that the zmax is sorted by value, so the highest peak is
        %first and the corresponding index is also first in imax, the same pattern
        %hold for zmin and imin
        sorted_max_indexes = sort(imax);
        highest_max_index = find(sorted_max_indexes == imax(1));
        
        %locate the index between the highest and next maximums
        min_index = find(imin > sorted_max_indexes(highest_max_index) & imin < sorted_max_indexes(highest_max_index + 1));
        assert(length(min_index) == 1, 'Error: expected to only find one minimum index');
        mask_thresh = intensity(imin(min_index));
    end
end

if (isnan(mask_thresh))
    mask_image_orig = double(imread(I_file));
    mask_min_max = csvread(fullfile(fileparts(I_file),filenames.raw_mask_min_max));
    mask_image_orig_norm = (mask_image_orig - min(mask_min_max))/range(mask_min_max);
    
    [pathstr,name, ~] = fileparts(I_file);
    out_file = fullfile(pathstr,[name,'_edge.png']);

    imwrite(mask_image_orig_norm,out_file);
    
    return
end

threshed_mask = mask_image > mask_thresh;

%%Mask Cleanup
connected_areas = bwlabel(threshed_mask);
region_sizes = regionprops(connected_areas, 'Area'); %#ok<MRPBW>

%filter out connected regions smaller than 2.5% of the total image area
threshed_mask = ismember(connected_areas, find([region_sizes.Area] > 0.025*length(mask_image(:))));

threshed_mask = imfill(threshed_mask,'holes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pixel Intensity Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hist(mask_pixels(:),length(unique(mask_pixels(:)))/4);
xlabel('Pixel Intensity');
ylabel('Pixel Count');
hold on;
ylimits = ylim;
plot([mask_thresh,mask_thresh],ylimits,'g','LineWidth',3);

if (exist('smoothed_heights','var'))
    plot(intensity,smoothed_heights,'r','LineWidth',3);
    plot(intensity(imax),zmax,'x');
end

hold off;
print('-dpng',fullfile(fileparts(I_file),'cell_mask_hist.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_file = fullfile(fileparts(I_file),filenames.cell_mask);

imwrite(threshed_mask, out_file);

if (i_p.Results.debug)
    addpath('../visualize_cell_features/');
    
    mask_image_orig = double(imread(I_file));
    mask_min_max = csvread(fullfile(fileparts(I_file),filenames.raw_mask_min_max));
    mask_image_orig_norm = (mask_image_orig - min(mask_min_max))/range(mask_min_max);
    
    edge_highlight = create_highlighted_image(mask_image_orig_norm,bwperim(threshed_mask),'color_map',[1,0,0]);
    [pathstr,name, ~] = fileparts(I_file);
    out_file = fullfile(pathstr,[name,'_edge.png']);
    
    imwrite(edge_highlight,out_file);
end