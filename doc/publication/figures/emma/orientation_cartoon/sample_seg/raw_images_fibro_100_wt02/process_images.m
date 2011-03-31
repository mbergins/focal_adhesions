adhesions = imread('adhesions_perim.png');
focal_image = double(imread('focal_image.png'));
focal_norm = (focal_image - min(focal_image(:)))/range(focal_image(:));

col_1 = 145;
col_2 = 382;
row_1 = 500;
row_2 = 699;

imwrite(focal_norm(row_1:row_2,col_1:col_2),'focal_norm.png')

ad_highlights = create_highlighted_image(focal_norm,adhesions > 0,'color_map',[0,0/255,255/255]);
imwrite(ad_highlights(row_1:row_2,col_1:col_2,:),'ad_highlights.png')

major_axis = csvread('raw_data/MajorAxisLength.csv');
minor_axis = csvread('raw_data/MinorAxisLength.csv');

eccen = major_axis./minor_axis;

high_eccen = eccen >= 3;
low_eccen = eccen < 3;

high_eccen_adhesions = ismember(adhesions, find(high_eccen));
low_eccen_adhesions = ismember(adhesions, find(low_eccen));

addpath('/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')

highlights = create_highlighted_image(focal_norm,high_eccen_adhesions,'color_map',[0,1,0]);
highlights = create_highlighted_image(highlights,low_eccen_adhesions,'color_map',[1,0,0]);

imwrite(highlights(row_1:row_2,col_1:col_2,:),'eccen_highlight.png')

only_high_eccen = create_highlighted_image(focal_norm,high_eccen_adhesions,'color_map',[0,1,0]);

imwrite(only_high_eccen(row_1:row_2,col_1:col_2,:),'high_eccen_highlight.png')