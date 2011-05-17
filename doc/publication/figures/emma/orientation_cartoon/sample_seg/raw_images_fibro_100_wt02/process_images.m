%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File Reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')

adhesions = imread('adhesions_perim.png');
adhesions_full = imread('adhesions.png');
focal_image = double(imread('focal_image.png'));
focal_norm = (focal_image - min(focal_image(:)))/range(focal_image(:));

col_1 = 145;
col_2 = 382;
row_1 = 500;
row_2 = 699;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Highlight Creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imwrite(focal_norm(row_1:row_2,col_1:col_2),'focal_norm.png')

ad_highlights = create_highlighted_image(focal_norm,adhesions > 0,'color_map',[105/255,249/255,255/255]);
imwrite(ad_highlights(row_1:row_2,col_1:col_2,:),'ad_highlights.png')
% 
% major_axis = csvread('raw_data/MajorAxisLength.csv');
% minor_axis = csvread('raw_data/MinorAxisLength.csv');
% orientation = csvread('raw_data/Orientation.csv');
% centroid = csvread('raw_data/Centroid.csv');
% 
% eccen = major_axis./minor_axis;
% 
% high_eccen = eccen >= 3;
% low_eccen = eccen < 3;
% 
% high_eccen_adhesions = ismember(adhesions, find(high_eccen));
% low_eccen_adhesions = ismember(adhesions, find(low_eccen));
% 
% addpath('/home/mbergins/Documents/Projects/focal_adhesions/trunk/src/visualize_cell_features/')
% 
% highlights = create_highlighted_image(focal_norm,high_eccen_adhesions,'color_map',[0,1,0]);
% highlights = create_highlighted_image(highlights,low_eccen_adhesions,'color_map',[1,0,0]);
% 
% imwrite(highlights(row_1:row_2,col_1:col_2,:),'eccen_highlight.png')
% 
% only_high_eccen = create_highlighted_image(focal_norm,high_eccen_adhesions,'color_map',[0,1,0]);
% imwrite(only_high_eccen(row_1:row_2,col_1:col_2,:),'high_eccen_highlight.png')
% 
% output_file = 'high_eccen_highlight_nc.png';
% imwrite(only_high_eccen,'high_eccen_highlight_nc.png')
% 
% high_ad_nums = unique(adhesions(high_eccen_adhesions));
% for i=1:length(high_ad_nums)
%     obj_num = high_ad_nums(i);
%     
%     pos_str = [' +',num2str(centroid(obj_num,1)),'+',num2str(centroid(obj_num,2))];
%     label_str = [' "',num2str(obj_num),'\n', ...
%         sprintf('%.0f',orientation(obj_num)),'"'];
%     command_str = ['convert ', output_file, ' -fill ''rgba(255,255,255)'' -annotate', ...
%         pos_str, label_str, ' ', output_file];
%     system(command_str);
% end
% 
% Pull out adhesion 172
% 
% ad_172 = adhesions == 172;
% ad_172_high = create_highlighted_image(focal_norm,ad_172,'color_map',[0,1,0],'mix_percent',0.5);
% limits = [543,563,330,340];
% imwrite(ad_172_high(limits(1):limits(2),limits(3):limits(4),:),'ad_172_high.png')
% 
% %Pull out adhesion 100
% ad_100 = adhesions == 100;
% ad_100_high = create_highlighted_image(focal_norm,ad_100,'color_map',[0,1,0],'mix_percent',0.5);
% limits = [650,672,276,287];
% 
% imwrite(ad_100_high(limits(1):limits(2),limits(3):limits(4),:),'ad_100_high.png')