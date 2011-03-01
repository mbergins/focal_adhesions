function data_set = read_in_file_set(this_dir,filenames)

data_set = struct;

data_set.this_dir = this_dir;

data_set.focal_image = imread(fullfile(this_dir,filenames.focal_image));