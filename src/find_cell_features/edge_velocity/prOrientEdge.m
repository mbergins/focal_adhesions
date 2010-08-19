function [pixel_list]= prOrinetEdge(img_mask, ORIENT)
% PRORIENTEDGE 
% 
%            
%
% SYNOPSIS        [pixel_list]=prOrientEdge(img_mask, ORIENT)
%
% INPUT           img_mask   :
%                 ORIENT     :
% 
% OUTPUT          pixel_list : oriented edge pixel list
%                           
% DEPENDENCES     prOrinetEdge uses (                                   
%                                   }
%                 prOrinetEdge is used by {imEdgeTracker }
%
% Matthias Machacek 22/03/06

%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[tmp_boundary] = bwboundaries(img_mask, 'noholes');
pixel_list(:,1) = tmp_boundary{1}(:,2);
pixel_list(:,2) = tmp_boundary{1}(:,1);
clear tmp_boundary;


if ORIENT == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start edge at the leftmost part of the edge              %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [dum, ind] = min(pixel_list(:,1));
    
elseif ORIENT == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start edge at the righttmost part of the edge           %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [dum, ind] = max(pixel_list(:,1));
elseif  ORIENT == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start edge at the highest part of the edge              %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    [dum, ind] = min(pixel_list(:,2));   
elseif  ORIENT == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start edge at the lowest part of the edge               %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [dum, ind] = max(pixel_list(:,2));      
elseif ORIENT == 10   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % orient pixels acording to enclosing ellipse orientation  %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [L,num] = bwlabel(img_mask,8);
    if num == 1
        cell_prop = regionprops(L,'Centroid','Orientation','MajorAxisLength','MinorAxisLength');
        cell_centroid = cell_prop.Centroid;
        cell_orientation = cell_prop.Orientation;
        cell_major_axis_l = cell_prop.MajorAxisLength;
        cell_minor_axis_l = cell_prop.MinorAxisLength;

        % minor axis position at ellipse boundary
        x_e = cell_centroid(1) + cell_minor_axis_l/2 * cos(90-cell_orientation);
        y_e = cell_centroid(2) - cell_minor_axis_l/2 * sin(90-cell_orientation);
        % get the closest pixel to the minor axis position

        % get the closest boundary pixel
        d = (pixel_list(:,1) - x_e).^2 + (pixel_list(:,2) - y_e).^2;
        [val, ind] = min(d);
        %figure,imshow(mask);
        %hold on
        %plot(cell_centroid(1),cell_centroid(2),'.r');
        %plot(x_e,y_e,'.r');
        %plot(pixel_list(ind,1),pixel_list(ind,2),'.');
    end
elseif  ORIENT >= 4 && ORIENT <= 7  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % orient pixels acording to a cross-hair set to the center of gravity %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [L,num] = bwlabel(img_mask,8);
    if num == 1
        cell_prop = regionprops(L,'Centroid','Orientation','MajorAxisLength','MinorAxisLength');
        cell_centroid = round(cell_prop.Centroid);
        if ORIENT == 4 
            ind_1 = find(cell_centroid(2) == pixel_list(:,2));
            [dum ind_2] = min(pixel_list(ind_1,1));
            ind = ind_1(ind_2);
        end
        if ORIENT == 5 
            ind_1 = find(cell_centroid(2) == pixel_list(:,2));
            [dum ind_2] = max(pixel_list(ind_1,1));
            ind = ind_1(ind_2);
        end    
        if ORIENT == 6
            ind_1 = find(cell_centroid(1) == pixel_list(:,1));
            [dum ind_2] = max(pixel_list(ind_1,2));
            ind = ind_1(ind_2);
        end
        if ORIENT == 7
            ind_1 = find(cell_centroid(1) == pixel_list(:,1));
            [dum ind_2] = min(pixel_list(ind_1,2));
            ind = ind_1(ind_2);
        end
    end
end %ORIENT


% figure
% plot(pixel_list(:,1),pixel_list(:,2),'.');
% hold on
% plot(cell_centroid(1),cell_centroid(2),'xr');
% plot(pixel_list(ind,1),pixel_list(ind,2),'xr');


% ind is new pixel number one, re-sort the others
temp_pixel_list = pixel_list;
clear pixel_list;
pixel_list_1 = temp_pixel_list(ind:end,:);
pixel_list_2 = temp_pixel_list(1:ind-1,:);
pixel_list = cat(1,pixel_list_1,pixel_list_2);
clear pixel_list_1;
clear pixel_list_2;
%figure
%plot(pixel_list(:,1),pixel_list(:,2));
    




