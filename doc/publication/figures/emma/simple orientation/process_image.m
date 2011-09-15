%% Image comes from Fibro_100ug_trial_2/WT_10 - Image 001

addpath('~/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/');

focal_image = double(imread('focal_image.png'));
adhesions = imread('adhesions.png');

limits = [62,469,70,564];

focal_image = (focal_image - min(focal_image(:)))/range(focal_image(:));

props = regionprops(adhesions,'MajorAxisLength','MinorAxisLength');

ratio = [props.MajorAxisLength]./[props.MinorAxisLength];

high_ratio = ismember(adhesions,find(ratio >= 3));
low_ratio = ismember(adhesions,find(ratio < 3));

highlight = create_highlighted_image(focal_image,high_ratio,'color_map',[0,1,0]);
highlight = create_highlighted_image(highlight,low_ratio,'color_map',[1,0,0]);

highlight = highlight(limits(1):limits(2),limits(3):limits(4),1:3);

imwrite(highlight,'eccen_above_3.png');