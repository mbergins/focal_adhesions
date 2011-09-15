addpath('~/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NS Stable Adhesion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('variable_angle_fibro_100ug_WT_01_ad1.mat');

ad_1 = all_images{1};

%remove final that follows death
ad_1 = ad_1(1:(length(ad_1)-1));

write_montage_image_set(ad_1,'stable_angle_full_set.png');

%trim set down to first frame and then a few intermediate
ad_1 = ad_1(1:6:end);
addpath('/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')

montage_ad_no_border = write_montage_image_set(ad_1,'stable_angle.png','num_cols',length(ad_1), ... 
    'pixel_size',0.1333333333,'bar_size',5);
imwrite(montage_ad_no_border,'NS_stable_angle.png')

b_size = 10;
montage_ad = write_montage_image_set(ad_1,'stable_angle.png','num_cols',length(ad_1), ... 
    'pixel_size',0.1333333333,'bar_size',5,'border_size',b_size,'border_color',[0,1,0]);

alphas = ones(size(montage_ad,1),size(montage_ad,2));
alphas(1:b_size,:) = 0.5;
alphas((end-b_size+1):end,:) = 0.5;
alphas(:,1:b_size) = 0.5;
alphas(:,(end-b_size+1):end) = 0.5;

imwrite(montage_ad,'sample_NS_stable_angle.png','Alpha',alphas)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2xKD Stable Adhesion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('variable_angle_100ug_pax13_ad1202.mat');

ad_1202 = temp;

while (isempty(ad_1202{1})), ad_1202 = ad_1202(2:end); end

%remove final that follows death
ad_1202 = ad_1202(2:end);
ad_1202 = ad_1202(1:(length(ad_1202)-1));

write_montage_image_set(ad_1202,'sample_2xKD_ad.png');

%trim set down to first frame and then a few intermediate
ad_1202 = ad_1202(1:6:end);
addpath('/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')

montage_ad_no_border = write_montage_image_set(ad_1202,'2xKD_set.png','num_cols',length(ad_1202), ... 
    'pixel_size',0.1333333333,'bar_size',5);
imwrite(montage_ad_no_border,'2xKD_stable_angle.png')

b_size = 10;
montage_ad = write_montage_image_set(ad_1202,'2xKD_set.png','num_cols',length(ad_1202), ... 
    'pixel_size',0.1333333333,'bar_size',5,'border_size',10,'border_color',[1,0,0]);

alphas = ones(size(montage_ad,1),size(montage_ad,2));
alphas(1:b_size,:) = 0.5;
alphas((end-b_size+1):end,:) = 0.5;
alphas(:,1:b_size) = 0.5;
alphas(:,(end-b_size+1):end) = 0.5;

imwrite(montage_ad,'sample_2xKD_stable_angle.png','Alpha',alphas)