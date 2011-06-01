function data_set = read_in_file_set(this_dir,filenames)

data_set = struct;

data_set.this_dir = this_dir;

data_set.focal_image = imread(fullfile(this_dir,filenames.focal_image));
data_set.focal_image_range = csvread(fullfile(this_dir,filenames.focal_image_min_max));

data_set.focal_norm = (double(data_set.focal_image) - data_set.focal_image_range(1))/ ... 
    range(data_set.focal_image_range);

data_set.adhesions = imread(fullfile(this_dir,filenames.adhesions));
data_set.adhesions_perim = imread(fullfile(this_dir,filenames.adhesions_perim));