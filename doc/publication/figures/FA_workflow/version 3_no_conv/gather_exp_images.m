focal_image = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/001/focal_image.png');
cell_mask = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/001/cell_mask.png');
ad_binary = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/001/adhesions_binary.png');
ad_label = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/001/adhesions.png');
temp_props = regionprops(ad_label,'area');

addpath('../../../../src/visualize_cell_features/')

big_edge = bwperim(imdilate(cell_mask,strel('disk',4)));

both = or(big_edge,ad_binary);

bounds = find_bounding_box(both);
bounds(1:2) = bounds(1:2) - 2;
bounds(3:4) = bounds(3:4) + 2;

small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));
small_both = small_both(1:200,200:end);

small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));
small_focal = small_focal(1:200,200:end);

imwrite(small_focal,'small_focal_1.png');
imwrite(not(small_both),'edge_and_ad_1.bmp');
system('/sw/bin/potrace -s -t 0 edge_and_ad_1.bmp');
system('rm edge_and_ad_1.bmp');

focal_image = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/012/focal_image.png');
cell_mask = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/012/cell_mask.png');
ad_binary = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/012/adhesions_binary.png');

addpath('../../../../src/visualize_cell_features/')

big_edge = bwperim(imdilate(cell_mask,strel('disk',4)));

both = or(big_edge,ad_binary);

small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));
small_both = small_both(1:200,200:end);

small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));
small_focal = small_focal(1:200,200:end);

imwrite(small_focal,'small_focal_2.png');
imwrite(not(small_both),'edge_and_ad_2.bmp');
system('/sw/bin/potrace -s -t 0 edge_and_ad_2.bmp');
system('rm edge_and_ad_2.bmp');

focal_image = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/024/focal_image.png');
cell_mask = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/024/cell_mask.png');
ad_binary = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/024/adhesions_binary.png');

addpath('../../../../src/visualize_cell_features/')

big_edge = bwperim(imdilate(cell_mask,strel('disk',4)));

both = or(big_edge,ad_binary);

small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));
small_both = small_both(1:200,200:end);

small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));
small_focal = small_focal(1:200,200:end);

imwrite(small_focal,'small_focal_3.png');
imwrite(not(small_both),'edge_and_ad_3.bmp');
system('/sw/bin/potrace -s -t 0 edge_and_ad_3.bmp');
system('rm edge_and_ad_3.bmp');

%%Second Try

% first_i_num = '001';
% focal_image = imread(fullfile('../../../../results/focal_adhesions/time_series_01/individual_pictures',first_i_num,'focal_image.png'));
% cell_mask = imread(fullfile('../../../../results/focal_adhesions/time_series_01/individual_pictures',first_i_num,'cell_mask.png'));
% ad_binary = imread(fullfile('../../../../results/focal_adhesions/time_series_01/individual_pictures',first_i_num,'adhesions_binary.png'));
% ad_label = imread(fullfile('../../../../results/focal_adhesions/time_series_01/individual_pictures',first_i_num,'adhesions.png'));
% 
% addpath('../../../../src/visualize_cell_features/')
% 
% big_edge = bwperim(imdilate(cell_mask,strel('disk',4)));
% 
% both = or(big_edge,ad_binary);
% 
% small_both = both(325:475,175:325);
% 
% small_focal = focal_image(325:475,175:325);
% 
% imwrite(small_focal,'small_focal_1.png');
% imwrite(not(small_both),'edge_and_ad_1.bmp');
% system('/sw/bin/potrace -s -t 0 edge_and_ad_1.bmp');
% system('rm edge_and_ad_1.bmp');
% 
% focal_image = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/026/focal_image.png');
% cell_mask = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/026/cell_mask.png');
% ad_binary = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/026/adhesions_binary.png');
% 
% addpath('../../../../src/visualize_cell_features/')
% 
% big_edge = bwperim(imdilate(cell_mask,strel('disk',4)));
% 
% both = or(big_edge,ad_binary);
% 
% small_both = both(325:475,175:325);
% 
% small_focal = focal_image(325:475,175:325);
% 
% imwrite(small_focal,'small_focal_2.png');
% imwrite(not(small_both),'edge_and_ad_2.bmp');
% system('/sw/bin/potrace -s -t 0 edge_and_ad_2.bmp');
% system('rm edge_and_ad_2.bmp');
% 
% focal_image = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/045/focal_image.png');
% cell_mask = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/045/cell_mask.png');
% ad_binary = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/045/adhesions_binary.png');
% 
% addpath('../../../../src/visualize_cell_features/')
% 
% big_edge = bwperim(imdilate(cell_mask,strel('disk',4)));
% 
% small_both = both(325:475,175:325);
% 
% small_focal = focal_image(325:475,175:325);
% 
% imwrite(small_focal,'small_focal_3.png');
% imwrite(not(small_both),'edge_and_ad_3.bmp');
% system('/sw/bin/potrace -s -t 0 edge_and_ad_3.bmp');
% system('rm edge_and_ad_3.bmp');


%%Third Try

% addpath('../../../../src/visualize_cell_features/')
% 
% first_i_num = '01';
% focal_image = imread(fullfile('../../../../results/focal_adhesions/time_series_04/individual_pictures',first_i_num,'focal_image.png'));
% cell_mask = imread(fullfile('../../../../results/focal_adhesions/time_series_04/individual_pictures',first_i_num,'cell_mask.png'));
% ad_binary = imread(fullfile('../../../../results/focal_adhesions/time_series_04/individual_pictures',first_i_num,'adhesions_binary.png'));
% ad_label = imread(fullfile('../../../../results/focal_adhesions/time_series_04/individual_pictures',first_i_num,'adhesions.png'));
% 
% big_edge = bwperim(imdilate(cell_mask,strel('disk',3)));
% 
% both = or(big_edge,ad_binary);
% 
% bounds = find_bounding_box(both);
% bounds(1:2) = bounds(1:2) - 2;
% bounds(3:4) = bounds(3:4) + 2;
% 
% small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));
% 
% small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));
% 
% imwrite(small_focal,'small_focal_1.png');
% imwrite(not(small_both),'edge_and_ad_1.bmp');
% system('/sw/bin/potrace -s -t 0 edge_and_ad_1.bmp');
% system('rm edge_and_ad_1.bmp');
% 
% focal_image = imread('../../../../results/focal_adhesions/time_series_04/individual_pictures/10/focal_image.png');
% cell_mask = imread('../../../../results/focal_adhesions/time_series_04/individual_pictures/10/cell_mask.png');
% ad_binary = imread('../../../../results/focal_adhesions/time_series_04/individual_pictures/10/adhesions_binary.png');
% 
% big_edge = bwperim(imdilate(cell_mask,strel('disk',3)));
% 
% both = or(big_edge,ad_binary);
% 
% small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));
% 
% small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));
% 
% imwrite(small_focal,'small_focal_2.png');
% imwrite(not(small_both),'edge_and_ad_2.bmp');
% system('/sw/bin/potrace -s -t 0 edge_and_ad_2.bmp');
% system('rm edge_and_ad_2.bmp');
% 
% focal_image = imread('../../../../results/focal_adhesions/time_series_04/individual_pictures/20/focal_image.png');
% cell_mask = imread('../../../../results/focal_adhesions/time_series_04/individual_pictures/20/cell_mask.png');
% ad_binary = imread('../../../../results/focal_adhesions/time_series_04/individual_pictures/20/adhesions_binary.png');
% 
% big_edge = bwperim(imdilate(cell_mask,strel('disk',3)));
% 
% both = or(big_edge,ad_binary);
% 
% small_both = both(bounds(2):bounds(4), bounds(1):bounds(3));
% 
% small_focal = focal_image(bounds(2):bounds(4), bounds(1):bounds(3));
% 
% imwrite(small_focal,'small_focal_3.png');
% imwrite(not(small_both),'edge_and_ad_3.bmp');
% system('/sw/bin/potrace -s -t 0 edge_and_ad_3.bmp');
% system('rm edge_and_ad_3.bmp');