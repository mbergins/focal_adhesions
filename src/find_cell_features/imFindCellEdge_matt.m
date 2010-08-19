function [ans, img_edge, img_bg, edge_pixel, length_edge, cell_isolated, rem_pix, object_position_old, varargout]=imFindCellEdge(img_org, img_name, contr, varargin)
% IMFINDCELLEDGE finds the the boundary of a object in an image 
% 
%             This function finds the leading edge of a cell (or any object in
%             in general) based on a pixel intensity criterion thus the
%             image is segmented by thresholding. The function 
%             "imFindThreshFilt" finds the correct threshold level. After
%             thresholding there might be several detected objects.
%             Thereafter the edge pixels of the object are extracted and
%             sorted.
%             The image intensities must be normalized to the range [0..1].
%             The image should also already be Gauss-filtered.
%             Currently not used:
%             Attention: if the image is filtered first (filter_image=1) it
%             will be also cropped to remove the filter artefacts on the image
%             boarders!
%            
%
% SYNOPSIS        [img_edge, img_bg, edge_pixel, length_edge, frame_pos]=imFindCellEdge(img_org, img_name, contr, varargin)
%
% INPUT           img_org:      image 
%                 img_name:     string identifier of the image
%                 contr:        flag for the display of control images
% 
% OUTPUT          ans:          status of success: 1: ok, -1:not ok
%                 img_edge:     image with the extracted edge
%                 img_bg:       mask: object: 1, background: 0
%                 edge_pixel:   2D matrix with the edge pixel
%                 length_edge:  the length of the border
%                 cell_isolated:  1 if cell is isolated 0 touches img border
%                 rem_pix:      pixels that were removed to eliminate noise
%                 object_position_old: I *think* this is the x,y coordinate
%                 of the center of the largest object in the image, but it
%                 appears not to be used; always returns 1
%                 varargs:
%                   P0:           the probability of the individual normal distributions
%                   MU0:          the mean of the individual normal distributions
%                           
% DEPENDENCES     imFindCellEdge uses ( Gauss2D,
%                                       imFindThreshFilt
%                                       bw_col_edge_filter,                                  
%                                   }
%                 imFindCellEdge is used by {imEdgeTracker }
%
% Matthias Machacek 25/09/03

%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin)
	l=length(varargin);
    for i=1:2:l
        in_found=0;
        if strcmp(varargin(i),'median_f')
            MEDIAN_F=varargin{i+1};  
            in_found=1;
        elseif strcmp(varargin(i),'normalize')
            NORMALIZE=varargin{i+1};  
            in_found=1;
        elseif strcmp(varargin(i),'bit_depth')
            BIT_DEPTH=varargin{i+1};  
            in_found=1;
        elseif strcmp(varargin{i},'edge_ pixel_nr')
            EDGE_PIXEL_NR=varargin{i+1}; 
            in_found=1;
        elseif strcmp(varargin{i},'img_sigma')
            IMG_SIGMA=varargin{i+1};     
            in_found=1;
        elseif strcmp(varargin{i},'filter_image')
            FILTER_IMAGE=varargin{i+1};    
            in_found=1;
        elseif strcmp(varargin{i},'erode_dilate')
            ERODE_DILATE=varargin{i+1};   
            in_found=1;
        %%% cluster parameters %%%%%%%%%%%%%%%%%%%    
        elseif strcmp(varargin{i},'cluster')
            CLUSTER=varargin{i+1};   
            in_found=1;     
        elseif strcmp(varargin{i},'cluster_method')
            CLUSTER_METHOD=varargin{i+1};   
            in_found=1;   
        elseif strcmp(varargin{i},'k_cluster')
            K_CLUSTER=varargin{i+1};   
            in_found=1; 
        elseif strcmp(varargin{i},'k_max')
            K_MAX=varargin{i+1};   
            in_found=1;
        elseif strcmp(varargin{i},'k_min')
            K_MIN=varargin{i+1};   
            in_found=1;
        elseif strcmp(varargin{i},'binning')
            BINNING=varargin{i+1};   
            in_found=1;
        elseif strcmp(varargin(i),'p0')
            P0=varargin{i+1};  
            in_found=1;
        elseif strcmp(varargin(i),'mu0')
            MU0=varargin{i+1};
            in_found=1;
        %%% threshold function parameters %%%%%%%%%%
        elseif strcmp(varargin{i},'f_window')
            F_WINDOW=varargin{i+1}; 
            in_found=1;
        elseif strcmp(varargin{i},'f_sigma')
            F_SIGMA=varargin{i+1};    
            in_found=1;
        elseif strcmp(varargin(i),'manual_level')  
            MANUAL_LEVEL = varargin{i+1};   
            in_found=1;  
        elseif strcmp(varargin(i),'cell_mode')  
            CELL_MODE= varargin{i+1};   
            in_found=1;               
        end
        
        if in_found == 0
            error_string = char(varargin(i));
            error(['Unknown input:   ' , error_string]);
        end 
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Default parameter of current function %%%%%%%%%%%%%%%%%%%%
if ~exist('NORMALIZE','var')
   %depth for the normalization of the image
   NORMALIZE=0;
end
if ~exist('BIT_DEPTH','var')
   %bit depth
   BIT_DEPTH=16;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Cluster specific parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('CLUSTER','var')
   %use cluster algorithm
   CLUSTER=0;
end
if ~exist('CLUSTER_METHOD','var')
   %method to use: 'kmeans' or 'em'
   CLUSTER_METHOD='kmeans';
end
if ~exist('K_CLUSTER','var')
   %number of clusters for kmeans method
   K_CLUSTER=3;
end
if ~exist('K_MAX','var')
   %maximum number of cluster for EM algorithm
   K_MAX=3;
end
if ~exist('K_MIN','var')
   %minimum number of cluster for EM algorithm
   K_MIN=3;
end
if ~exist('CELL_MODE','var')
   % number of modes that are counted as cell
   CELL_MODE = 1;
end
if ~exist('BINNING','var')
   %binning of the image 
   BINNING=0;
end
%the probability of the individual normal distributions
if ~exist('P0','var')
   P0=[];
end
%the mean values of the individual normal distributions
if ~exist('MU0','var')
   MU0=[];
end
if ~exist('MANUAL_LEVEL','var')
    % manual level, if -1 it is not used
    MANUAL_LEVEL = -1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('MEDIAN_F','var')
   %box width of the median filter
   MEDIAN_F=3;
end
if ~exist('EDGE_PIXEL_NR','var')
   %dito with the edge !! it is not used anymore!!
   EDGE_PIXEL_NR=500;
end
if ~exist('IMG_SIGMA','var')
   %the sigma of the gauss filter used for
   %image smoothing
   IMG_SIGMA=0.9;
end
if ~exist('FILTER_IMAGE','var')
   %the sigma of the gauss filter used for
   %image smoothing
   FILTER_IMAGE=1;
end
if ~exist('ERODE_DILATE','var')
   %Errode Dilate the binary image in order 
   %to remove noise on the boundary
   ERODE_DILATE=6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters for the thresholding function "imFindThreshFilt"  %%
if ~exist('F_WINDOW','var')
   F_WINDOW=5;%11
end
if ~exist('F_SIGMA','var')
   F_SIGMA=0.1;%0.6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert raw mask to binary mask
% This is from find_cell_mask.m and find_adhesion_properties.m

% find the cell edge
img_edge = bwperim(im2bw(img_org,adaptive_thresh(img_org,'upper_mean_weight',0.2)));
[img_bg, img_edge] = clean_up_mask_image(img_edge);
%edge_pixel = ind2sub(size(img_edge),find(bwperim(img_edge)));

edge_labels = bwlabel(img_edge, 8);
edge_labels_stat = regionprops(edge_labels, 'PixelList');
edge_pixel = edge_labels_stat(1).PixelList;

% This is from imFindCellEdge.m
%if ERODE_DILATE
    %this is done to REMOVE NOISE on the edge
    %and to close open holes by "missing" speckles
    %close the image
%    img_bg_er_dil = imclose(img_bg,strel('disk',ERODE_DILATE));  
       
    %control how many pixels where removed
%    diff_img=img_bg-img_bg_er_dil;
%    rem_pix = abs(sum(sum(diff_img)));

%    img_bg=img_bg_er_dil;
%else
    rem_pix=0;
%end

% get the length of the edge
% This is from imFindCellEdge.m
edge_nr = size(edge_pixel);
sq = sqrt(2);
length_edge = 0;
for i = 2:edge_nr
    if edge_pixel(i,1) ~= edge_pixel(i-1,1) & edge_pixel(i,2) ~= edge_pixel(i-1,2)
        length_edge=length_edge + sq;
    else
        length_edge=length_edge + 1;
    end
end

% test if cell is isolated or touches the image border
% This is from imFindCellEdge.m
[n_img, m_img]=size(img_org);
bi1=find(img_edge(1:n_img,2), 1);
bi2=find(img_edge(1:n_img,m_img-1), 1);
bi3=find(img_edge(2,1:m_img), 1);
bi4=find(img_edge(n_img-1,1:m_img), 1);

if isempty(bi1) & isempty(bi2) & isempty(bi3) & isempty(bi4)
    cell_isolated = 1;
else
    cell_isolated = 0;        
end

object_position_old = 1;

varargout{1} = P0;   
varargout{2} = MU0;

ans = 1

end