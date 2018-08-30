function filenames = add_filenames_to_struct(filenames)
%ADD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                           various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenames.individual_results_dir = 'individual_pictures';

filenames.focal_image = 'focal_image.png';
filenames.focal_image_min_max = '../../adhesion_props/image_analysis/focal_min_max.csv';
filenames.focal_image_threshold = '../../adhesion_props/image_analysis/focal_filter_results.csv';

filenames.adhesions = 'adhesions.png';
filenames.adhesions_binary = 'adhesions_binary.png';
filenames.adhesions_perim = 'adhesions_perim.png';

filenames.raw_mask = 'raw_mask_file.png';
filenames.cell_mask = 'cell_mask.png';
filenames.raw_mask_min_max = '../../adhesion_props/image_analysis/cell_mask_min_max.csv';

filenames.focal_image_secondary = 'focal_image_secondary.png';
filenames.focal_image_secondary_min_max = '../../adhesion_props/image_analysis/focal_secondary_min_max.csv';
filenames.focal_image_secondary_threshold = '../../adhesion_props/image_analysis/focal_secondary_filter_results.csv';

filenames.kinase = 'kinase.png';
filenames.kinase_min_max = '../../adhesion_props/image_analysis/kinase_min_max.csv';

filenames.photo_bleach_correction = 'photo_correct.txt';

filenames.background_intensity = '../../adhesion_props/image_analysis/background_intensity.txt';

filenames.focal_image_threshold_plot = '../../adhesion_props/image_analysis/focal_threshold.eps';
filenames.per_image_threshold_plot = '../../adhesion_props/image_analysis/per_image_threshold.eps';

filenames.assembly_rows = '../../adhesion_props/assem_disassem_models/assembly_rows_lengths.csv';
filenames.disassembly_rows = '../../adhesion_props/assem_disassem_models/disassembly_rows_lengths.csv';

filenames.centroid_x = '../../adhesion_props/lin_time_series/Centroid_x.csv';
filenames.centroid_y = '../../adhesion_props/lin_time_series/Centroid_y.csv';

filenames.max_intensity = '../../visualizations/max_intensity_projection.png';

filenames.tracking = 'tracking_matrices/tracking_seq.csv';