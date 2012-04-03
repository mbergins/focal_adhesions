function filenames = add_filenames_to_struct(filenames)
%ADD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                           various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenames.focal_image = 'focal_image.png';

filenames.adhesions = 'adhesions.png';
filenames.adhesions_binary = 'adhesions_binary.png';
filenames.adhesions_perim = 'adhesions_perim.png';

filenames.raw_mask = 'raw_mask_file.png';
filenames.cell_mask = 'cell_mask.png';

filenames.photo_bleach_correction = 'photo_correct.txt';
filenames.background_intensity = '../../adhesion_props/background_intensity.txt';

filenames.focal_image_min_max = '../../adhesion_props/focal_min_max.csv';
filenames.focal_image_threshold = '../../adhesion_props/focal_filter_results.csv';
filenames.focal_image_threshold_plot = '../../adhesion_props/focal_threshold.eps';
filenames.per_image_threshold_plot = '../../adhesion_props/per_image_threshold.eps';

filenames.assembly_rows = '../../adhesion_props/assembly_rows_lengths.csv';
filenames.disassembly_rows = '../../adhesion_props/disassembly_rows_lengths.csv';

filenames.centroid_x = '../../adhesion_props/lin_time_series/Centroid_x.csv';
filenames.centroid_y = '../../adhesion_props/lin_time_series/Centroid_y.csv';

filenames.max_intensity = '../../visualizations/max_intensity_projection.png';
