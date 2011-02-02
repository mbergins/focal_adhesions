labeled_ad = imread('../../../../../results/focal_adhesions/time_series_01/individual_pictures/001/adhesions.png');
focal_image = double(imread('../../../../../results/focal_adhesions/time_series_01/individual_pictures/001/focal_image.png'));
focal_image = focal_image/2^16;

small_labeled_ad = labeled_ad(100:250,125:300);
small_focal_image = focal_image(100:250,125:300);

%renumber the found adhesions to start at one
ad_nums = unique(small_labeled_ad);
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(ad_nums)
    small_labeled_ad(small_labeled_ad == ad_nums(i)) = i - 1;
end

addpath(genpath('../../../../../src/visualize_cell_features/'))

c_map = jet(length(unique(small_labeled_ad)) - 1);
c_map = c_map(randperm(size(c_map,1)),:);

highlighted_ad = create_highlighted_image(small_focal_image,small_labeled_ad,'color_map', c_map);
imwrite(highlighted_ad,'pixel_ID/highlighted_ad.png')
imwrite(small_focal_image,'pixel_ID/small_sample_focal.png')