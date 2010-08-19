addpath('../../../../src/visualize_cell_features/')

base_focal_image = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/001/focal_image.png');
base_focal_image = double(base_focal_image)/2^16;

adhesions = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/001/adhesions_perim.png');
b_box = find_binary_bounding_box(adhesions);
b_box(1:2) = b_box(1:2) - 10;
b_box(3:4) = b_box(3:4) + 10;

cmap = jet(double(max(adhesions(:))));
cmap = cmap(randperm(size(cmap,1)),:);

highlighted_ad = create_highlighted_image(base_focal_image,adhesions,'color_map',cmap);

bounded_fi = base_focal_image(b_box(2):b_box(4), b_box(1):b_box(3),:);
bounded_ha = highlighted_ad(b_box(2):b_box(4), b_box(1):b_box(3),:);

bounded_fi = draw_scale_bar(bounded_fi,0.215051);
% bounded_ha = draw_scale_bar(bounded_ha,0.215051);

imwrite(bounded_fi,'focal_image.png')

imwrite(bounded_ha,'highlighted.png')