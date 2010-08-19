function [ans, img_edge, img_bg, edge_pixel, length_edge, frame_pos, rem_pix, varargout]=imFindCellEdge(img_org, img_name, contr, varargin)
% IMFINDCELLEDGE finds the the boundary of a object in an image 
% 
%             This function finds the leading edge of a cell (or any object in
%             in general) based on a pixel intensity criterion thus the
%             image is segmented by thresholding. The function 
%             "imFindThreshFilt" finds the correct threshold level. After
%             thresholding there might be several detected objects. There
%             right one is assumed to be the larger than (OBJ_SIZE).
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
%                 length_edge:  the length of the boarder
%                 frame_pos :   position of edge at image border
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
        elseif strcmp(varargin(i),'obj_size')
            OBJ_SIZE=varargin{i+1};
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
if ~exist('OBJ_SIZE','var')
   %the minimal area of a object in pixel
   %if object is smaller it will not be detected
   OBJ_SIZE=8000; %it used to be 25000
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

%check if the image is normalized
if (max(max(img_org)) > 1 | min(min(img_org)) < 0) &  NORMALIZE==0
    error('Image is not normalized to [0..1]');
end
if NORMALIZE
    if NORMALIZE == 8 | NORMALIZE == 10 | NORMALIZE == 12 | NORMALIZE == 14 | NORMALIZE == 16 | NORMALIZE == 255 ...
        | NORMALIZE == 1023 | NORMALIZE == 4095 | NORMALIZE == 16383 | NORMALIZE == 65635 
        
        if NORMALIZE < 24
            img_org=double(img_org)./(2^NORMALIZE-1);
            MANUAL_LEVEL = double(MANUAL_LEVEL)./(2^NORMALIZE-1);
        else 
            img_org = double(img_org)./NORMALIZE;
            MANUAL_LEVEL = double(MANUAL_LEVEL)./NORMALIZE;
        end
    else 
        error('unsuported image depth, only images of 8,10,12,14,16 bit are accepted.');
    end
end


%determine size of image
[n_img, m_img]=size(img_org);

if FILTER_IMAGE
    %gauss filtering of the image 
    %img_org=Gauss2D(img_org,IMG_SIGMA);

    img_org=Gauss2DBorder(img_org,IMG_SIGMA);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Threshold image  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CLUSTER
    if ~isempty(MU0)
        %use initial values
        [seg_img, P0, MU0]=imClusterSeg(img_org, contr, 'method', CLUSTER_METHOD,... 
                                        'k_cluster', K_CLUSTER, 'k_min', K_MIN, 'k_max', K_MAX, 'binning', BINNING,...
                                        'p0', P0, 'mu0', MU0, 'emptyaction', 'singleton');  %em
    else    
        %use random seeding
        [seg_img, P0, MU0]=imClusterSeg(img_org, contr, 'method', CLUSTER_METHOD,... 
                                        'k_cluster', K_CLUSTER, 'k_min', K_MIN, 'k_max', K_MAX, 'binning', BINNING,...
                                        'emptyaction', 'singleton');  %em                                  
    end
    ans = 1;
    img_thresh = seg_img == 1;

    %give the cluster poistions so they can be used in the next time step
    varargout{1} = P0;   
    varargout{2} = MU0;
    
    %obj_val is acctually not needed!
 
elseif MANUAL_LEVEL >= 0
        img_thresh = img_org <= MANUAL_LEVEL;
        thresh = MANUAL_LEVEL;
        ans = 1;
else
    if 0
        [minThresh, J] = imMinimumThreshold (img_org, BIT_DEPTH);
        thresh = minThresh;
        ans = 1;
    else
        [ans, thresh, int_bg, int_obj]=imFindThreshFilt(img_org, BIT_DEPTH, contr, 'f_window', F_WINDOW, 'f_sigma', F_SIGMA);
        %save the threshold for the next timestep
    end
    img_thresh=img_org < thresh;
    %here the object is BLACK (=0) and the background is WHITE (=1)!!!
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   End threshold image  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ans == -1
    %threshold failed by first method. Use clustering method!
    CLUSTER = 1;
    if ~isempty(MU0)
        display('Normal threshold failed, try cluster method with initial value');
        [ans, img_edge, img_bg, edge_pixel, length_edge, frame_pos, rem_pix, P0, MU0]=...
            imFindCellEdge(img_org, img_name, contr, 'bit_depth',BIT_DEPTH, 'filter_image',FILTER_IMAGE,  'img_sigma',IMG_SIGMA,...
            'f_window', F_WINDOW,  'f_sigma', F_SIGMA,...
            'erode_dilate',ERODE_DILATE, 'median_f', MEDIAN_F,...
            'cluster', CLUSTER, 'cluster_method', CLUSTER_METHOD, 'k_cluster', K_CLUSTER, 'k_min', K_MIN,...
            'k_max', K_MAX, 'p0', P0, 'mu0', MU0);
    else
        display('Normal threshold failed, try cluster method without initial value');
        [ans, img_edge, img_bg, edge_pixel, length_edge, frame_pos, rem_pix, P0, MU0]=...
            imFindCellEdge(img_org, img_name, contr, 'bit_depth',BIT_DEPTH, 'filter_image',FILTER_IMAGE,  'img_sigma',IMG_SIGMA,...
            'f_window', F_WINDOW,  'f_sigma', F_SIGMA,...
            'erode_dilate',ERODE_DILATE, 'median_f',MEDIAN_F,...
            'cluster', CLUSTER, 'cluster_method', CLUSTER_METHOD, 'k_cluster', K_CLUSTER, 'k_min', K_MIN,...
            'k_max', K_MAX);                                                                            
    end  
    %give the cluster positions so they can be used in the next time step
    varargout{1} = P0;   
    varargout{2} = MU0;
    return;
else
    %threshold succeded!
	if contr
        figure
        title(img_name);   
        subplot(2,2,1);
        imshow(img_thresh,[]);
        title('Thresholded image'); 
        if ~CLUSTER
            text(10,100,['Treshold value: ',num2str(thresh)],'Color','r'); 
        end
	end
	
	%do a median filtering
    warning off MATLAB:conversionToLogical;
	img=medfilt2(img_thresh,[MEDIAN_F MEDIAN_F]);
	
	%find all objects in the image
	img=~img;
	img_labels = bwlabel(img,8);
	labels_stat=regionprops(img_labels,'Area','PixelList');
	img_area=[labels_stat.Area];
	
    %get the largest object
    [img_area_rev, j] = max(img_area);
    
    %img_area_rev=img_area > OBJ_SIZE;
	%%find index of relevant region
	%[i,j,revelant]=find(img_area_rev);
    
	nr_found_obj=size(j,2);
	if nr_found_obj<1
        disp('No object found, check OBJ_SIZE');
        disp('But we try it with the clustering!');
        %if no solution was found use the clustering
        CLUSTER = 1;
        if ~isempty(MU0)
            display('No object found after threshold, try cluster method with initial value');
            [ans, img_edge, img_bg, edge_pixel, length_edge, frame_pos, rem_pix, P0, MU0]=...
                imFindCellEdge(img_org, img_name, contr, 'bit_depth',BIT_DEPTH, 'filter_image',FILTER_IMAGE,  'img_sigma',IMG_SIGMA,...
                'f_window', F_WINDOW,  'f_sigma', F_SIGMA,...
                'erode_dilate',ERODE_DILATE, 'median_f', MEDIAN_F,...
                'cluster', CLUSTER, 'cluster_method', CLUSTER_METHOD, 'k_cluster', K_CLUSTER, 'k_min', K_MIN,...
                'k_max', K_MAX, 'p0', P0, 'mu0', MU0);
        else
            display('No object found after threshold, try cluster method without initial value');
            [ans, img_edge, img_bg, edge_pixel, length_edge, frame_pos, rem_pix, P0, MU0]=...
                imFindCellEdge(img_org, img_name, contr, 'bit_depth',BIT_DEPTH, 'filter_image',FILTER_IMAGE,  'img_sigma',IMG_SIGMA,...
                'f_window', F_WINDOW,  'f_sigma', F_SIGMA,...
                'erode_dilate',ERODE_DILATE, 'median_f',MEDIAN_F,...
                'cluster', CLUSTER, 'cluster_method', CLUSTER_METHOD, 'k_cluster', K_CLUSTER, 'k_min', K_MIN,...
                'k_max', K_MAX);                                                                            
        end  
        %give the cluster poistions so they can be used in the next time step
        varargout{1} = P0;   
        varargout{2} = MU0;
        return;
	end
	img_ref(1:n_img,1:m_img)=j(1);
	img_bg=img_labels==img_ref;
	%fill the holes in the object (e.g. cell), background has value 0
	img_bg=imfill(double(img_bg),'holes');
	
	if ERODE_DILATE
        %this is done to REMOVE NOISE on the edge
        %and to close open holes by "missing" speckles
        %close the image
        img_bg_er_dil = imclose(img_bg,strel('disk',ERODE_DILATE));  
       
        %control how many pixels where removed
        diff_img=img_bg-img_bg_er_dil;
        rem_pix = abs(sum(sum(diff_img)));

        img_bg=img_bg_er_dil;
	else
        rem_pix=0;
	end
	
	%ok now you have to fill pseudo holes, thus holes
	%connected to the background
	img_labels_bg = bwlabel(~logical(img_bg),8);
	labels_stat_bg=regionprops(img_labels_bg,'Area');
	img_bg_area=[labels_stat_bg.Area];
    
    %test if there was any object found
    if isempty(img_bg_area)
        %if no solution was found use the clustering
        CLUSTER = 1;
        if ~isempty(MU0)
            display('No object found after threshold, try cluster method with initial value');
            [ans, img_edge, img_bg, edge_pixel, length_edge, frame_pos, rem_pix, P0, MU0]=...
                imFindCellEdge(img_org, img_name, contr, 'bit_depth',BIT_DEPTH, 'filter_image',FILTER_IMAGE,  'img_sigma',IMG_SIGMA,...
                'f_window', F_WINDOW,  'f_sigma', F_SIGMA,...
                'erode_dilate',ERODE_DILATE, 'median_f', MEDIAN_F,...
                'cluster', CLUSTER, 'cluster_method', CLUSTER_METHOD, 'k_cluster', K_CLUSTER, 'k_min', K_MIN,...
                'k_max', K_MAX, 'p0', P0, 'mu0', MU0);
        else
            display('No object found after threshold, try cluster method without initial value');
            [ans, img_edge, img_bg, edge_pixel, length_edge, frame_pos, rem_pix, P0, MU0]=...
                imFindCellEdge(img_org, img_name, contr, 'bit_depth',BIT_DEPTH, 'filter_image',FILTER_IMAGE,  'img_sigma',IMG_SIGMA,...
                'f_window', F_WINDOW,  'f_sigma', F_SIGMA,...
                'erode_dilate',ERODE_DILATE, 'median_f',MEDIAN_F,...
                'cluster', CLUSTER, 'cluster_method', CLUSTER_METHOD, 'k_cluster', K_CLUSTER, 'k_min', K_MIN,...
                'k_max', K_MAX);                                                                            
        end  
        %give the cluster poistions so they can be used in the next time step
        varargout{1} = P0;   
        varargout{2} = MU0;
        return;
    end
    
    
    [val_bg bg_index]=max(img_bg_area);
	img_bg=bg_index==img_labels_bg;
	
	if contr
        subplot(2,2,2);
        imshow(img_bg,[]);
        title('Median filtered, object size filtered'); 
        text(10,150, ['Median filter size: ',num2str(MEDIAN_F)],'Color','r');   
        text(10,100,['Obj. size threshold: ',num2str(OBJ_SIZE)],'Color','r');    
        if ERODE_DILATE
            subplot(2,2,3);
            imshow(~diff_img);
            title('Original-Eroded/Dilated image');
            text(10,100,['Number of removed pix: ',num2str(rem_pix)],'Color','r');        
        end
	end
	
    %get the edge
    img_edge = bwperim(img_bg,4);
    %clean the boarder
    img_edge(1:n_img,1)=0;
    img_edge(1:n_img,m_img)=0;
    img_edge(1,1:m_img)=0;
    img_edge(n_img,1:m_img)=0;

    % eliminate false detected edges based on a size criterion
    edge_labels2=bwlabel(img_edge,8);
    labels_edge_stat2=regionprops(edge_labels2,'Area','PixelList');
    % take just the largest edge
    edge_area_rev2=[labels_edge_stat2.Area] == max([labels_edge_stat2.Area]);
    % edge_area_rev2=[labels_edge_stat2.Area] > EDGE_PIXEL_NR;
    % find index (=label) of relevant edge
    [i2,j2,revelant2]=find(edge_area_rev2);
    if length(j2) > 1
        error('more than one object edge found, check EDGE_PIXEL_NR');
    elseif  length(j2)==0
        error('no object edge found, check EDGE_PIXEL_NR');
    end
    img_edge=edge_labels2==j2(1);

    %get the pixels
    edge_pixel=labels_edge_stat2(j2(1)).PixelList;

    
    % sort the pixels
    edge_pixel = imSortPixels(edge_pixel);


    % get the length of the edge
    [edge_nr dim]=size(edge_pixel);
    sq=sqrt(2);
    length_edge=0;
    for i=2:edge_nr
        if edge_pixel(i,1) ~= edge_pixel(i-1,1) & edge_pixel(i,2)~= edge_pixel(i-1,2)
            length_edge=length_edge+sq;
        else
            length_edge=length_edge+1;
        end
    end
    
    % get the position of the edge at the image border
    bi1=find(img_edge(1:n_img,2)==1);
    if size(bi1,2)>1 bi1=bi1'; end
    bi2=find(img_edge(1:n_img,m_img-1)==1);
    if size(bi2,2)>1 bi2=bi2'; end
    bi3=find(img_edge(2,1:m_img)==1);
    if size(bi3,2)>1 bi3=bi3'; end
    bi4=find(img_edge(n_img-1,1:m_img)==1);
    if size(bi4,2)>1 bi4=bi4'; end
    frame_pos = [bi1', bi2', bi3', bi4'];
    
    % test if object does not touch a boarder
    if isempty(frame_pos)
        % set variable to "undefined" status
        frame_pos=[0, 0];
    elseif length(frame_pos) > 2
        dum = frame_pos;
        clear frame_pos
        frame_pos = dum(1:2);
    end

		
	if contr
        %superimpose the original image and the extracted edge
        subplot(2,2,4);
        imshow(img_org,[]);
        hold on
        plot(edge_pixel(:,1), edge_pixel(:,2),'r');   
        title('Original image with extracted edge');
	end
	
	ans=1;
    img_bg=~img_bg;
    
    if ~CLUSTER
        %these are not needed but they are in the output
        varargout{1}=[];
        varargout{2}=[];
    end
end