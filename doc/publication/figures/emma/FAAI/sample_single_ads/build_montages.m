%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable Adhesion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('variable_angle_fibro_100ug_WT_01_ad506.mat');

ad_506 = all_images{506};

%remove empties and pre-birth image
ad_506 = ad_506(8:end);

%remove final that follows death
ad_506 = ad_506(1:(length(ad_506)-1));

write_montage_image_set(ad_506,'variable_angle_full_set.png');

%trim set down to first frame and then a few intermediate
ad_506 = ad_506(1:6:end);

addpath('/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')

write_montage_image_set(ad_506,'variable_angle.png','num_cols',length(ad_506), ... 
    'pixel_size',0.1333333333,'bar_size',5);

clear all_images;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stable Adhesion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('variable_angle_fibro_100ug_WT_01_ad1.mat');

ad_1 = all_images{1};

%remove final that follows death
ad_1 = ad_1(1:(length(ad_1)-1));

write_montage_image_set(ad_1,'stable_angle_full_set.png');

%trim set down to first frame and then a few intermediate
ad_1 = ad_1(1:6:end);
addpath('/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')

write_montage_image_set(ad_1,'stable_angle.png','num_cols',length(ad_1), ... 
    'pixel_size',0.1333333333,'bar_size',5);