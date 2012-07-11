function varargout = find_cell_mask(I_file,varargin)
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

if (not(any(strcmp(i_p.UsingDefaults,'mask_threshold'))))
    threshed_mask = mask_image > i_p.Results.mask_threshold;
    varargout{1} = i_p.Results.mask_threshold;
else
    %%Threshold identification
    
    [heights, intensity] = hist(mask_image(:),1000);
    
    smoothed_heights = smooth(heights,0.05,'loess');
    [~,imax,~,imin]= extrema(smoothed_heights);
    
    %keep in mind that the zmax is sorted by value, so the highest peak is
    %first and the corresponding index is also first in imax, the same pattern
    %hold for zmin and imin
    
    sorted_max_indexes = sort(imax);
    first_max_index = find(sorted_max_indexes == imax(1));
    
    %locate the index between the first two maximums
    min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
    assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));
    threshed_mask = mask_image > intensity(imin(min_index));
    
    hist(mask_image(:),1000);
    xlim([0,10000]);
    ylim([0,14000]);
    xlabel('Pixel Intensity');
    ylabel('Pixel Count');
    hold on;
    plot(intensity,smoothed_heights,'r','LineWidth',3);
    hold on;
    ylimits = ylim;
    plot([intensity(imin(min_index)),intensity(imin(min_index))],ylimits,'g','LineWidth',3);
    hold off;
    print('-dpng',fullfile(fileparts(I_file),'cell_mask_hist.png'))
    varargout{1} = intensity(imin(min_index));
end

%%Mask Cleanup
connected_areas = bwlabel(threshed_mask);
region_sizes = regionprops(connected_areas, 'Area');

%filter out connected regions smaller than 2.5% of the total image area
threshed_mask = ismember(connected_areas, find([region_sizes.Area] > 0.025*length(mask_image(:))));

threshed_mask = imfill(threshed_mask,'holes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_file = fullfile(fileparts(I_file),filenames.cell_mask);

imwrite(threshed_mask, out_file)