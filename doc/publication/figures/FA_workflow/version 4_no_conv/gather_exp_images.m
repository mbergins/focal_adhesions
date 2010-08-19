addpath('../../../../src/visualize_cell_features/')

base_dir = '../../../../results/focal_adhesions/time_series_01/individual_pictures';

focal_image = imread(fullfile(base_dir,'012/focal_image.png'));
cell_mask = imread(fullfile(base_dir,'012/cell_mask.png'));
ad_binary = imread(fullfile(base_dir,'012/adhesions_binary.png'));
ad_label = imread(fullfile(base_dir,'012/adhesions.png'));

big_edge = bwperim(imdilate(cell_mask,strel('disk',3)));

both = or(big_edge,ad_binary);

bounds = find_bounding_box(both);
bounds(1:2) = bounds(1:2) - 2;
bounds(3:4) = bounds(3:4) + 2;
bounds = [100,200,250,350];

small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));

small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));

small_focal = draw_scale_bar(small_focal,0.215051);

imwrite(small_focal,'small_focal_1.png');
imwrite(not(small_both),'edge_and_ad_1.bmp');
system('/sw/bin/potrace -s -t 0 edge_and_ad_1.bmp');
system('rm edge_and_ad_1.bmp');

focal_image = imread(fullfile(base_dir,'030/focal_image.png'));
cell_mask = imread(fullfile(base_dir,'030/cell_mask.png'));
ad_binary = imread(fullfile(base_dir,'030/adhesions_binary.png'));
ad_label = imread(fullfile(base_dir,'030/adhesions.png'));

big_edge = bwperim(imdilate(cell_mask,strel('disk',3)));

both = or(big_edge,ad_binary);

small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));

small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));

imwrite(small_focal,'small_focal_2.png');
imwrite(not(small_both),'edge_and_ad_2.bmp');
system('/sw/bin/potrace -s -t 0 edge_and_ad_2.bmp');
system('rm edge_and_ad_2.bmp');

focal_image = imread(fullfile(base_dir,'048/focal_image.png'));
cell_mask = imread(fullfile(base_dir,'048/cell_mask.png'));
ad_binary = imread(fullfile(base_dir,'048/adhesions_binary.png'));
ad_label = imread(fullfile(base_dir,'048/adhesions.png'));

big_edge = bwperim(imdilate(cell_mask,strel('disk',3)));

both = or(big_edge,ad_binary);

small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));

small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));

imwrite(small_focal,'small_focal_3.png');
imwrite(not(small_both),'edge_and_ad_3.bmp');
system('/sw/bin/potrace -s -t 0 edge_and_ad_3.bmp');
system('rm edge_and_ad_3.bmp');