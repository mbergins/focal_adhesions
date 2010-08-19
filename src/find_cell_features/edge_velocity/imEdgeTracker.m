function [img_proccessed, img_edge]=imEdgeTracker(varargin)
% IMEDGETRACKER measures the protrusion rate of cell leading edge
%
%               imEdgeTracker measures the protrusion rate of a cell edge. For
%               that a time series of cell images, such as FSM, is
%               required. The image is segmented based on an intensity
%               criterion. The threshold value for this is automaticly determined
%               based on the image histogram.
%               The segmented binary image is median filtered and
%               cleaned up. Attention, if the cell covers only a small part
%               of the image unexpected results can occure since a size
%               criterion is used to identify the cell body.
%               The cell edge is defined as 8-nghb. connected. For the
%               analysis the edge is approximated with a smoothing
%               spline.
%               The edge protrusion (displacement) is calculated.
%               There are various definitions of the protrusion, i)
%               displacement along the normal to the leading edge, ii)
%               along the shortest distance between two positions and iii)
%               based on a mechanical model.
%               The mechanical model needs a good estimate as initial value,
%               preferabely from the "nearest" calcualtion!
%               After stack processing some data analyis is performed.
%               As a control output the variable 'img_proccessed' is given.
%               The values are: 0 if the image is not processed (because
%               t_step is not 1), 1 for a successful process and -1 when
%               the atempt failed
%               If you want to calculate the protrusion:
%               'protrusion',1
%
%               For the protrusion vector calculation you can define a edge
%               indent from the left and right edge side. Do this because
%               the edge might be detected badly at the edges and this
%               might disturb the protrusion calculation. However the
%               segments start at the edge start and end at the edge end.
%               Therefore the beginning and ending segments might contain
%               no protrusion vectors. In this case the average vectors get
%               the value -99
%
%               Parameters
%               'bit_depth',16383
%               'max_img',20
%               'contr',1
%               'movie',0
%               'write_data',0
%               'lambda',500
%               'na',1.4
%               'bit_depth',16383
%               'max_img',5
%               't_step',1
%               'protrusion',1
%               'prot_sampling',10
%               'prot_depth',20
%               'nr_sect',10
%               'filter_image',1
%               'img_sigma',0.9
%               'normal',0
%               'nearest',1
%               'mechanical',1
%               'k_s',0.1
%               'k_w',1
%
%               The function writes the following data to disk:
%               the coordinates of the edge pixels
%               object mask (binary image)
%               averaged unit normal vectors
%               averaged protrusion vectors
%               edge_spline.mat:    the spline time series cell edge description
%
% SYNOPSIS      [img_edge]=imEdgeTracker()
%
% INPUT
%
% OUTPUT                   img_proccessed : -1 error, 0 not processed, 1 ok
%                          img_edge       : gives the image with the edge
%
%
%
% DEPENDENCES   imEdgeTracker uses { uigetfile,
%                                    imreadstacknd,
%                                    getFileStackNames,
%                                    imFindCellEdge,
%                                    imFindThreshFilt,
%                                    imPixelChainSpline,
%                                    prGetDisplNormal,
%                                    prGetDispNearest,
%                                    prGetDispMech,
%                                    prGetProtRegion,
%                                    prGetAvEdge,
%                                   }
%
% Matthias Machacek 10/14/03

%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);
for i=1:2:l
    in_found=0;
    if strcmp(varargin(i),'contr')
        CONTR=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'debug')
        DEBUG=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'file')
        FILE=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'results')
        RESULTS=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'movie')
        MOVIE=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'pixel')
        PIXEL=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'time_interval')
        TIME_INTERVAL=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'lambda')
        LAMBDA=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'na')
        NA=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'bit_depth')
        BIT_DEPTH=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'first_img')
        FIRST_IMG=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'max_img')
        MAX_IMG=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'t_step')
        T_STEP=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'protrusion')
        PROTRUSION=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'prot_sampling')
        PROT_SAMPLING=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'prot_depth')
        PROT_DEPTH=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'nr_sect')
        NR_SECT=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'filter_image')
        FILTER_IMAGE=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'img_sigma')
        IMG_SIGMA=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'filter_f')
        FILTER_F=varargin{i+1};
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(varargin{i},'f_window')
        F_WINDOW=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin{i},'f_sigma')
        F_SIGMA=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin{i},'tolerance')
        TOLERANCE=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'normal')
        NORMAL=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'nearest')
        NEAREST=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'robust_min')
        ROBUST_MIN=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'tol')
        TOL=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'mechanical')
        MECHANICAL=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'level_set')
        LEVEL_SET=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'k_S')
        K_S=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'k_W')
        K_W=varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'parenth_l')
        PARENTH_L = varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'parenth_r')
        PARENTH_R = varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'movie')
        MOVIE = varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'use_bw_mask')
        USE_BW_MASK = varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'manual_level')
        MANUAL_LEVEL = varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'orient')
        ORIENT_CELL = varargin{i+1};
        in_found=1;
    elseif strcmp(varargin(i),'cell_mode')
        CELL_MODE = varargin{i+1};
        in_found=1;
    end

    if in_found == 0
        error_string = char(varargin(i));
        error(['Unknown input:   ' , error_string]);
    end
end

if ~exist('CONTR','var')
    %parameter for displaying control images
    CONTR=1;
end
if ~exist('USE_BW_MASK','var')
    % parameter for using data from pre run segmentation
    USE_BW_MASK = 0;
end
if ~exist('DEBUG','var')
    %parameter for debuging
    DEBUG=0;
end
if ~exist('FILE','var')
    %path to image data
    FILE=0;
end
if ~exist('RESULTS','var')
    %path to folder for the results
    RESULTS=0;
end
if ~exist('FIRST_IMG','var')
    %path to image data
    FIRST_IMG=1;
end
if ~exist('MOVIE','var')
    %parameter for displaying control images
    MOVIE=0;
end
if ~exist('LAMBDA','var')
    %the wave lenght used for the microscop imaging
    LAMBDA=500;
end
if ~exist('NA','var')
    %the apperture used for the microscop imaging
    NA=1.4;
end
if ~exist('PIXEL','var')
    %pixel size in nm
    PIXEL=67;
end
if ~exist('TIME_INTERVAL','var')
    %time interval between the image frames [s]
    TIME_INTERVAL=10;
end
if ~exist('BIT_DEPTH','var')
    %the bit depth of the images
    BIT_DEPTH=16;
end
if ~exist('PROTRUSION','var')
    %flag to include the protrusion calculation
    PROTRUSION=0;
end
if ~exist('PROT_SAMPLING','var')
    %protrusion sampling definded as "every .. pixel"
    PROT_SAMPLING=5;
end
if ~exist('PROT_DEPTH','var')
    %protrusion band mask width
    PROT_DEPTH=20;
end
if ~exist('NR_SECT','var')
    %protrusion band sector number
    NR_SECT=20;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters for "imFindCellEdge"  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('FILTER_IMAGE','var')
    %(Gauss) filtering of the image
    FILTER_IMAGE=1;
end
if ~exist('IMG_SIGMA','var')
    % Gauss filter variance
    IMG_SIGMA=0.9;
end
if ~exist('MEDIAN_F','var')
    %width of the median filter
    MEDIAN_F=3;
end
if ~exist('ERODE_DILATE','var')
    %gives the radius of the structuring element
    ERODE_DILATE=6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters for the thresholding function "imFindThreshFilt"  %%
if ~exist('F_WINDOW','var')
    %window width of the filter for image histogram
    F_WINDOW=5;
end
if ~exist('F_SIGMA','var')
    %inverse!! variance of filter for image histogram
    F_SIGMA=0.1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('MANUAL_LEVEL','var')
    % manual level, if -1 it is not used
    MANUAL_LEVEL = -1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Cluster specific parameters  used in "imClusterSeg" %%%%%
if ~exist('CLUSTER','var')
    %use cluster algorithm
    CLUSTER=0;
end
if ~exist('CLUSTER_METHOD','var')
    %use cluster algorithm
    CLUSTER_METHOD='kmeans';
end
if ~exist('K_CLUSTER','var')
    %use cluster algorithm
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters for the spline function "imPixelChainSpline"  %%%%%%
if ~exist('TOLERANCE','var')
    %smoothing parameter of the approximating spline
    TOLERANCE=30;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters for "prGetDispMech" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('NORMAL','var')
    %protrusion definition
    NORMAL=0;
end
if ~exist('NEAREST','var')
    NEAREST=0;
end
if ~exist('MECHANICAL','var')
    MECHANICAL=1;
end
if ~exist('LEVEL_SET','var')
    LEVEL_SET = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters for "prGetDispNearest" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ROBUST_MIN','var')
    % flag for discrete robust method
    ROBUST_MIN=1;
end
if ~exist('TOL','var')
    % spline parameter search spsce. If zero then entire spline is searched
    TOL=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter of spring constants in "prGetDispMech" %%%%%%%%%%%%%%%
if ~exist('K_S')
    %spacing
    K_S=0.1;
end
if ~exist('K_W')
    %angle
    K_W=1;
end
if ~exist('PARENTH_L')
    % set the parenthesis from the image border
    PARENTH_L = 1;
end
if ~exist('PARENTH_R')
    % set the parenthesis from the image border
    PARENTH_R = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter of how to determine beginning of edge for whole cells%
if ~exist('ORIENT_CELL')
    % orient = 0: start left side
    % orient = 1: start right side
    % orient = 2: start upper side
    % orient = 3: start lower side
    % orient = 4: ellipse
    ORIENT_CELL = 1;
end

BIT_DEPTH = 2 ^ BIT_DEPTH -1;



%LEVELSET = 0;
PATHS = 0;
%MECHANICAL = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% End parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if a image/ result path was supplied. If not, ask the user
% for that
if FILE
    firstfilename=FILE;
else
    [fileName,dirName, dialog_ans] = uigetfile('*.tif','Select image');
    if dialog_ans == 0
        img_proccessed = 0;
        img_edge       = 0;
        return
    end
    firstfilename=[dirName,fileName];
end
if RESULTS
    dir_w=RESULTS;
else
    dir_w = uigetdir('start_path','Select/create directory for results')
    if dir_w == 0
        img_proccessed = 0;
        img_edge       = 0;
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% get alternative image for overlay   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alt_img_path = '/public/forMatthias/mg0144_488_crop';
filelist_tmp = dir([alt_img_path filesep '*.tif']);

iEntry = 1;
for i=1:length(filelist_tmp)
    if(~filelist_tmp(i).isdir)
        filelist_alt_img(iEntry,:) = filelist_tmp(i).name;
        iEntry = iEntry + 1;
    end
end
%%%%%%%% get alternative image for overlay   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% create the protrusion directory within the working directory
% suppress the existing directory warning
warning off all
exist_dir = exist(dir_w,'dir');
if exist_dir == 0
    %it does not exist, so try to create it
    [s, mess, messid] = mkdir(dir_w);
    if s==0
        %creation failed, return with error message
        img_proccessed = 0;
        img_edge       = 0;
        disp('Failed to create the protrusion directory');
        return
    end
end

% create all the sub directories
mkdir(dir_w, 'cell_mask');
% folder for the cell, edge overlays
mkdir(dir_w, 'edge_cell');
% folder for the edge - vectors
mkdir(dir_w, 'pr_vectors');
% folder for the result figures
mkdir(dir_w, 'figures');

% enable warnings again
warning on all

if ~exist('MAX_IMG','var')
    % the index of the last image. If t_step is equal to one this is the
    % max numerber of images
    ans = inputdlg('Number of images','Number of images');
    MAX_IMG = str2num(ans{1});
    clear ans;
end
if ~exist('T_STEP','var')
    % step size
    ans = inputdlg('Time step','Time step');
    T_STEP = str2num(ans{1});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Write all the parameters to a file  %%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([dir_w filesep 'parameters.dat'],'w+');
if fid == -1
    error('Could not create file in results directory. Be sure to put a "edge" directory into results directory!!');
end

fprintf(fid, 'CONTR = %d\n',            CONTR);
fprintf(fid, 'FILE = %s\n',             firstfilename);
fprintf(fid, 'RESULTS = %s\n',          dir_w);
fprintf(fid, 'FIRST_IMG = %d\n',        FIRST_IMG);
fprintf(fid, 'MAX_IMG = %d\n',          MAX_IMG);
fprintf(fid, 'T_STEP = %d\n',           T_STEP);
fprintf(fid, 'MOVIE = %d\n',            MOVIE);
fprintf(fid, 'BIT_DEPTH = %d\n',        BIT_DEPTH);
fprintf(fid, 'LAMBDA = %d\n',           LAMBDA);
fprintf(fid, 'NA = %d\n',               NA);
fprintf(fid, 'PIXEL = %d\n',            PIXEL);
fprintf(fid, 'TIME_INTERVAL = %d\n',    TIME_INTERVAL);
fprintf(fid, 'PROTRUSION = %d\n',       PROTRUSION);
fprintf(fid, 'PROT_SAMPLING = %d\n',    PROT_SAMPLING);
fprintf(fid, 'PROT_DEPTH = %d\n',       PROT_DEPTH);
fprintf(fid, 'NR_SECT = %d\n',          NR_SECT);
fprintf(fid, 'FILTER_IMAGE = %d\n',     FILTER_IMAGE);
fprintf(fid, 'IMG_SIGMA = %d\n',        IMG_SIGMA);
fprintf(fid, 'MEDIAN_F = %d\n',         MEDIAN_F);
fprintf(fid, 'F_WINDOW = %d\n',         F_WINDOW);
fprintf(fid, 'F_SIGMA = %d\n',          F_SIGMA);
fprintf(fid, 'ERODE_DILATE = %d\n',     ERODE_DILATE);
fprintf(fid, 'TOLERANCE = %d\n',        TOLERANCE);
fprintf(fid, 'NORMAL = %d\n',           NORMAL);
fprintf(fid, 'NEAREST = %d\n',          NEAREST);
fprintf(fid, 'MECHANICAL = %d\n',       MECHANICAL);
fprintf(fid, 'TOL = %d\n',              TOL);
fprintf(fid, 'ROBUST_MIN = %d\n',       ROBUST_MIN);
fprintf(fid, 'PARENTH_L = %d\n',        PARENTH_L);
fprintf(fid, 'PARENTH_R = %d\n',        PARENTH_R);
fprintf(fid, 'K_S = %d\n',              K_S);
fprintf(fid, 'K_W = %d\n',              K_W);
fprintf(fid, 'CLUSTER = %d\n',          CLUSTER);
fprintf(fid, 'CLUSTER_METHOD = %s\n',   CLUSTER_METHOD);
fprintf(fid, 'K_CLUSTER = %d\n',        K_CLUSTER);
fprintf(fid, 'K_MAX = %d\n',            K_MAX);
fprintf(fid, 'K_MIN = %d\n',            K_MIN);
fprintf(fid, 'LEVEL_SET = %d\n',        LEVEL_SET);
fprintf(fid, 'MANUAL LEVEL = %d\n',     MANUAL_LEVEL);

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open a log file, store any information here!
% [fid, full_log_name]=createLogFile(firstfilename,'numbering','automatic');
% fprintf(fid,'Number of images: %i\n',MAX_IMG);
% fprintf(fid,'Time step: %i  image(s)\n',T_STEP);
% fclose(fid);

if MOVIE
    MakeQTMovie('start', [dir_w 'prot_movie.mov']);
end

% create a colormap for the movie edges
cmap_edge_evolution=jet(MAX_IMG);


%%%%%%%%%%%%%%%%%%    Allocate memory  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine the image dimension
[n_img, m_img]=size(imread(firstfilename));
% get the filelist filenames
[filelist]=getFileStackNames(firstfilename);

disp('Allocating memory');
img_edge_rgb=zeros(n_img,m_img,3);
img_edge_w_rgb=ones(n_img,m_img,3);
img_edge_sup=zeros(n_img,m_img);
img_proccessed=zeros(MAX_IMG*T_STEP,1);
rem_pix=zeros(MAX_IMG,1);

if PROTRUSION
    PROTRUSION_ARRAY_SIZE = 700;
    %disp_mech = zeros(PROTRUSION_ARRAY_SIZE,2,ceil(MAX_IMG));
    %edge_mech = zeros(PROTRUSION_ARRAY_SIZE,2,ceil(MAX_IMG));
    %disp_pos  = zeros(PROTRUSION_ARRAY_SIZE,ceil(MAX_IMG));
end
%%%%%%%%%%%%%%%%%% End allocate memory  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Start Time Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extracting cell edges');
index=0;
bit_depth_test = 1;

strg = sprintf('%%.%dd',3);
backSpc = ['\b\b\b'];

for time=FIRST_IMG: T_STEP: FIRST_IMG + T_STEP*(MAX_IMG-1)
    index=index+1;

    fprintf(1,[strg],index);

    fileName=char(filelist(time));
    [tmp_path,tmp_fname] = fileparts(fileName);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Image Segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check the bit depth of the image
    info = imfinfo(fileName);
    if info.BitDepth ~= round(log(BIT_DEPTH+1)/log(2)) & bit_depth_test
        ans = questdlg('The acctual bit depth is different from the user specified one. Continue?','Warning');
        if strcmp(ans,'No')
            img_proccessed = 0;
            img_edge       = 0;
            return
        else
            bit_depth_test =0;
        end
    end

    % read the image and normalize the intensity to a [0..1] range
    img_org=imreadnd2(fileName,0,BIT_DEPTH);

    % subtract background approximation
    % [img_org, avgBackLevel] = imSubtractBackgroundRegression(img_org);

    % find the cell edge
    if isempty(MU0)
        P0 = [];
        MU0 = [];
    end
    [ans, img_edge(:,:), mask(:,:), pixel_list, edge_l(index),...
        cell_isolated(index), rem_pix(index), cell_pos_old, P0, MU0]=...
        imFindCellEdge(img_org,fileName,CONTR, 'bit_depth',BIT_DEPTH, 'filter_image',FILTER_IMAGE,  'img_sigma',IMG_SIGMA,...
        'f_window', F_WINDOW,  'f_sigma', F_SIGMA,...
        'erode_dilate',ERODE_DILATE, 'median_f', MEDIAN_F,...
        'cluster', CLUSTER, 'cluster_method', CLUSTER_METHOD, 'k_cluster', K_CLUSTER, 'k_min', K_MIN,...
        'k_max', K_MAX, 'p0', P0, 'mu0', MU0, 'manual_level', MANUAL_LEVEL, 'cell_mode', CELL_MODE);

    % In case, re-orient pixel_list
    if index > 1
        d1 = sqrt((pixel_list(1,1) - pixel_list_last(1,1))^2+(pixel_list(1,2) - pixel_list_last(1,2))^2);
        d2 = sqrt((pixel_list(1,1) - pixel_list_last(end,1))^2+(pixel_list(1,2) - pixel_list_last(end,2))^2);
        if d2 < d1
            pixel_list = flipdim(pixel_list,1);
        end
    end
    % if cell is isolated, re-sort the pixels
    if cell_isolated(index) == 1
        [pixel_list]= prOrientEdge(mask, ORIENT_CELL);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% End image Segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ans== -1
        disp('No edge extracted');
        img_proccessed(time)=-1;
    else
        img_proccessed(time)= 1;

        pixel_edge{index} = pixel_list;

        imwrite(mask,fullfile(dir_w, 'cell_mask',  ['mask_' tmp_fname '.tif']),'tif');


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Create an image with the edge derived from the spline   %%%%
        [x_dim y_dim]=size(pixel_list);
        % create a image of superimposed edges
        img_edge_sup=img_edge_sup+img_edge(:,:);
        % get the spline describing the cell edge
        % the parameter runs from [1..number of data points]
        [edge_sp_x, edge_sp_y]=imPixelChainSpline(pixel_list(:,:),'tolerance',TOLERANCE);
        edge_sp_array_x(index)=edge_sp_x;
        edge_sp_array_y(index)=edge_sp_y;

        % store the number of pixel of the edge
        edge_pix_nr(index)=max(size(pixel_list));


        % get the edge pixels described by the spline
        p_n=1:edge_pix_nr(index);
        spline_pixel=[fnval(edge_sp_x,p_n); fnval(edge_sp_y,p_n)]';

        %superimpose the original image with the extracted edge
        if 0 %read alternative image as defined in line 476
            img_overlay = imread([alt_img_path filesep filelist_alt_img(time,:)]);
        else
            img_overlay = img_org;
        end
        %I'm not sure what nrm does, so I'll comment it out and hope
        %nothing breaks - MEB 9/25/2008
        %         sup_img = nrm(img_overlay,1);
        sup_img = img_overlay;
        sup_img = cat(3,sup_img,sup_img,sup_img);


        if time ==  FIRST_IMG
            sup_img_edge_evolution = sup_img;
            sup_img_edge_evolution_spline = sup_img;
            h_sup_img_edge_evolution_spline = figure('Visible','Off');
            imshow(sup_img_edge_evolution_spline);
        end

        for i=1:x_dim
            sup_img(pixel_list(i,2),pixel_list(i,1),3)=1;
            bool_y = round(spline_pixel(i,2)) >= 1 & round(spline_pixel(i,2)) <= n_img;
            bool_x = round(spline_pixel(i,1)) >= 1 & round(spline_pixel(i,1)) <= m_img;
            if bool_y & bool_x
                % this is single frame overlay
                sup_img(round(spline_pixel(i,2)),round(spline_pixel(i,1)),1) = 1;
                sup_img(round(spline_pixel(i,2)),round(spline_pixel(i,1)),2) = 0;
                sup_img(round(spline_pixel(i,2)),round(spline_pixel(i,1)),3) = 0;
                % this is integrated overlay with color coded edges
                sup_img_edge_evolution(round(spline_pixel(i,2)),round(spline_pixel(i,1)),1) = cmap_edge_evolution(index,1);
                sup_img_edge_evolution(round(spline_pixel(i,2)),round(spline_pixel(i,1)),2) = cmap_edge_evolution(index,2);
                sup_img_edge_evolution(round(spline_pixel(i,2)),round(spline_pixel(i,1)),3) = cmap_edge_evolution(index,3);
            end
        end
        warning off all;
        imwrite(sup_img, fullfile(dir_w, 'edge_cell', ['img_edge_' tmp_fname '.tif']),'tif');
        warning on all;

        % this is integrated overlay with spline edges
        figure(h_sup_img_edge_evolution_spline);
        set(h_sup_img_edge_evolution_spline,'Visible','Off');
        hold on
        plot(spline_pixel(:,1),spline_pixel(:,2));
        %%% END create an image with the edge derived from the spline  %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Get normal directions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create an array of spline parameters
        clear p_n;
        p_n=1:3:edge_pix_nr(index);

        % extract a band behind the protrusion edge and the
        % unit normal pointing away from the cell
        [x_normal_out, y_normal_out]=prGetProtRegion(img_org, edge_sp_x, edge_sp_y, p_n, 'prot_depth', PROT_DEPTH);

        % average the unit normals
        [x_normal_out_av, y_normal_out_av, x_n_av, y_n_av, m_pos]=...
            prGetAvEdge(edge_sp_x, edge_sp_y, edge_pix_nr(index), p_n, x_normal_out, y_normal_out, 'nr_sect', NR_SECT);

        % re-nornamize the averaged normals
        l_n_av = sqrt(x_normal_out_av.^2 + y_normal_out_av.^2);
        x_normal_out_av = x_normal_out_av ./ l_n_av;
        y_normal_out_av = y_normal_out_av ./ l_n_av;


        x_n = fnval(edge_sp_x, p_n);
        y_n = fnval(edge_sp_y, p_n);

        normal_matrix{index} = [x_n', y_n', x_normal_out', y_normal_out'];
        %%%%%%%%% End get normal directions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Operations on images involving two time steps %%%%%%%%%%%%%
        if index > 1
            if PROTRUSION
                % Set the points at which the protrusion is to be calculated
                i_nn=1:PROT_SAMPLING:edge_last_pix_nr;
                l_i_nn=length(i_nn);
                i_nn(l_i_nn-PARENTH_R:l_i_nn)=[];
                i_nn(1:PARENTH_L)=[];
                clear l_i_nn;

                if NORMAL
                    [temp1, temp2, i_pos]=prGetDispNormal(edge_sp_x_last, edge_sp_y_last, edge_sp_x, edge_sp_y, i_nn);

                    nr_prot_vectors(index) = size(temp1,1);
                    if nr_prot_vectors(index) > PROTRUSION_ARRAY_SIZE
                        error('Too many protrusion vectors! Increase array size!')
                    end
                    edge_mech(1:nr_prot_vectors(index) ,: ,index)   = temp1;
                    disp_mech(1:nr_prot_vectors(index) ,: ,index)   = temp2;
                    disp_pos(1:nr_prot_vectors(index), index)       = i_pos';
                    clear temp1 temp2 i_pos;

                elseif NEAREST
                    [temp1, temp2, i_pos]=prGetDispNearest(edge_sp_x_last, edge_sp_y_last, edge_sp_x, edge_sp_y, i_nn,...
                        'tol', TOL, 'robust_min', ROBUST_MIN);

                    nr_prot_vectors(index) = size(temp1,1);
                    if nr_prot_vectors(index) > PROTRUSION_ARRAY_SIZE
                        error('Too many protrusion vectors! Increase array size!')
                    end
                    edge_mech(1:nr_prot_vectors(index) ,: ,index)   = temp1;
                    disp_mech(1:nr_prot_vectors(index) ,: ,index)   = temp2;
                    disp_pos(1:nr_prot_vectors(index), index)       = i_pos';
                    clear temp1 temp2 i_pos;

                elseif MECHANICAL
                    % first find a initial solution based on the nearest
                    % model
                    [temp1, temp2, i_n]=prGetDispNearest(edge_sp_x_last, edge_sp_y_last, edge_sp_x, edge_sp_y, i_nn,...
                        'tol', TOL, 'robust_min', ROBUST_MIN);
                    clear temp1 temp2;

                    % take just the part of the spline with valid
                    % displacements. This is because the object can move out
                    % of the image
                    i_n_l = length(i_n);
                    if i_n(i_n_l) > edge_pix_nr(index)
                        i = i_n_l;
                        while i_n(i) > edge_pix_nr(index) & length(i_n) > 0
                            i_n(i) = [];
                            i_nn(i) = [];
                            i = i-1;
                        end
                    end
                    if i_n(1) < 1
                        i = 1;
                        while i_n(i) < 1 & length(i_n) > 0
                            i_n(i) =  [];
                            i_nn(i) = [];
                        end
                    end

                    % now use it for the mechanical model
                    i_0 = i_n;
                    clear i_n;
                    % calculate the protrusion based on a mechanical model
                    if length(i_0) > 2
                        [temp1, temp2, i_pos, x_normal, y_normal]=...
                            prGetDispMechFixL(edge_sp_x_last, edge_sp_y_last, edge_sp_x, edge_sp_y, i_nn, i_0, CONTR,...
                            'k_S', K_S, 'k_W', K_W);
                        nr_prot_vectors(index) = size(temp1,1);

                        protrusion{index-1} = [temp1 temp2];

                        clear temp1 temp2 i_pos;
                    else
                        nr_prot_vectors(index) = 1;
                        edge_mech(1:nr_prot_vectors(index) ,: ,index) = 0;
                        disp_mech(1:nr_prot_vectors(index) ,: ,index) = 0;
                        disp_pos(1:nr_prot_vectors(index), index)     = 0;
                    end
                elseif LEVEL_SET
                    % use the Level Set method
                    if index == 2
                        % This is for the first time step, here we have to
                        % initialize both level set functions
                        % This function gives two outputs:
                        % disp_points   : corresponding point at next timestep
                        % protrusion    : integrated marker path length
                        [phi_t1, val_t1, disp_points, ls_protrusion] =...
                            lsLineMatching('mask_img_t0',mask_last, 'mask_img_t1',mask,...
                            'x_spline_t0',edge_sp_x_last, 'y_spline_t0',edge_sp_y_last,...
                            'x_spline_t1',edge_sp_x, 'y_spline_t1',edge_sp_y,...
                            'known_zero_level_points_t0', pixel_list_last,...
                            'known_zero_level_points_t1', pixel_list,...
                            'i_nn',i_nn,'result_dir',dir_w,...
                            'control',1);
                    else
                        [phi_t1, val_t1, disp_points, ls_protrusion] =...
                            lsLineMatching('mask_img_t0',mask_last, 'mask_img_t1',mask,...
                            'x_spline_t0',edge_sp_x_last, 'y_spline_t0',edge_sp_y_last,...
                            'x_spline_t1',edge_sp_x, 'y_spline_t1',edge_sp_y,...
                            'known_zero_level_points_t0', pixel_list_last,...
                            'known_zero_level_points_t1', pixel_list,...
                            'val_t0', val_t0, 'phi_t0', phi_t0,...
                            'i_nn',i_nn,'result_dir',dir_w,...
                            'control',1);
                    end

                    x_nn=fnval(edge_sp_x_last,i_nn);
                    y_nn=fnval(edge_sp_y_last,i_nn);

                    nr_prot_vectors(index) = size(ls_protrusion,1);
                    edge_mech(1:nr_prot_vectors(index) ,: ,index)= [x_nn', y_nn'];
                    if 1
                        % linear interpolation
                        disp_mech(1:nr_prot_vectors(index) ,: ,index)= [disp_points(:,1) - x_nn', disp_points(:,2) - y_nn'];
                    else
                        % integrated paths, that does not work as it is!
                        disp_mech(1:nr_prot_vectors(index) ,: ,index) = abs(ls_protrusion);
                    end
                    disp_pos(1:nr_prot_vectors(index), index) = i_nn';

                    clear disp_points;
                    clear ls_protrusionl
                end  %if NEAREST, NORMAL, MECHANICAL, LEVEL_SET

                % save control images
                if 1
                    h_prot_control = figure('Visible','Off');
                    quiver(protrusion{index-1}(:,1), protrusion{index-1}(:,2), protrusion{index-1}(:,3),protrusion{index-1}(:,4), 0);
                    hold on
                    plot(fnval(edge_sp_x,      1: edge_sp_x.knots(end)),        fnval(edge_sp_y,      1: edge_sp_y.knots(end)));
                    plot(fnval(edge_sp_x_last, 1: edge_sp_x_last.knots(end)),   fnval(edge_sp_y_last, 1: edge_sp_y_last.knots(end)));
                    axis equal
                    axis ij

                    axis([0 m_img 0 n_img]);
                    axis off
                    %set(gca,'xtick',[]);
                    %set(gca,'ytick',[]);

                    print(h_prot_control,  fullfile(dir_w, 'pr_vectors', ['prot_vec_' tmp_fname '.tif']),'-dtiff');
                end
            end %if PROTRUSION
        end %if index>1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % if ans== -1

    %%%%%%%%%%   save current timestep variables for next timestep %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %img_org_last=img_org;
    pixel_list_last=[];
    pixel_list_last=pixel_list;
    edge_last_pix_nr=[];
    edge_last_pix_nr = edge_pix_nr(index);
    edge_sp_x_last=[];
    edge_sp_y_last=[];
    edge_sp_x_last=edge_sp_x;
    edge_sp_y_last=edge_sp_y;
    mask_last = mask;
    if index >1 & LEVEL_SET
        val_t0 = val_t1;
        phi_t0 = phi_t1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %make a rgb image
    [p q]=size(pixel_list);
    for j=1:p
        img_edge_rgb(pixel_list(j,2),pixel_list(j,1),1)=cmap_edge_evolution(index,1);
        img_edge_rgb(pixel_list(j,2),pixel_list(j,1),2)=cmap_edge_evolution(index,2);
        img_edge_rgb(pixel_list(j,2),pixel_list(j,1),3)=cmap_edge_evolution(index,3);

        img_edge_w_rgb(pixel_list(j,2),pixel_list(j,1),1)=cmap_edge_evolution(index,1);
        img_edge_w_rgb(pixel_list(j,2),pixel_list(j,1),2)=cmap_edge_evolution(index,2);
        img_edge_w_rgb(pixel_list(j,2),pixel_list(j,1),3)=cmap_edge_evolution(index,3);
    end


    if MOVIE
        MakeQTMovie('addmatrix',img_edge_rgb);
    end
    % save the last image of the movie
    if time ==  FIRST_IMG + T_STEP*(MAX_IMG-1)
        fileName=char(filelist(FIRST_IMG)); %T_STEP*(MAX_IMG-1)
        img_org=imreadnd2(fileName,0,BIT_DEPTH);
        imwrite(img_edge_rgb, [dir_w  '/figures' filesep 'edges_overlay.tif'],'tif' );
        imwrite(img_edge_w_rgb, [dir_w  '/figures' filesep 'edges_w_overlay.tif'],'tif');
        h_img_edge_rgb = figure;
        set(h_img_edge_rgb,'Visible','Off');
        imshow(img_edge_rgb);
        print(h_img_edge_rgb, [dir_w  '/figures' filesep 'edges_overlay.eps'],'-depsc2','-tiff');
        close(h_img_edge_rgb);
        h_img_edge_w_rgb = figure;
        set(h_img_edge_w_rgb,'Visible','Off');
        imshow(img_edge_w_rgb);
        print(h_img_edge_w_rgb, [dir_w  '/figures' filesep 'edges_w_overlay.eps'],'-depsc2','-tiff');
        close(h_img_edge_w_rgb);
    end
    clear pixel_list;

    % reset counter in matlab display
    fprintf(1,backSpc);
end

disp(['Total number of processed images: ', num2str(length(img_proccessed))]);
disp(['Successful processed images:      ', num2str(sum(img_proccessed))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% save the spline time series of the edge %%%%%%
save(fullfile(dir_w, 'edge_spline.mat'), 'edge_sp_array_x', 'edge_sp_array_y');
save(fullfile(dir_w, 'pixel_edge.mat'), 'pixel_edge');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PROTRUSION
    save(fullfile(dir_w, 'protrusion.mat'),    'protrusion');
end
save(fullfile(dir_w, 'normal_matrix.mat'), 'normal_matrix');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Postprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot superposition of protrusion and edge %%%%%%%%%%%%
if PROTRUSION
    h_prot_edge = figure;
    imshow(img_edge_rgb);
    hold on
    for ii=2:index
        quiver(protrusion{index-1}(:,1), protrusion{index-1}(:,2), protrusion{index-1}(:,3), protrusion{index-1}(:,4), 0,'r');
    end
    title('Cell leading edge and protrusion vectors');
    text(100,100,['Number of images  :',num2str(index)],'Color','r');
    text(100,120,['Image increment   :',num2str(T_STEP)],'Color','r');
    text(100,140,['Images  :',fileName],'Color','r','Interpreter','none');

    h_prot_edge_w = figure;
    imshow(img_edge_w_rgb);
    hold on
    for ii=2:index
        quiver(protrusion{index-1}(:,1), protrusion{index-1}(:,2), protrusion{index-1}(:,3), protrusion{index-1}(:,4), 0,'r');
    end
    title('Cell leading edge and protrusion vectors');
    text(100,100,['Number of images  :',num2str(index)],'Color','r');
    text(100,120,['Image increment   :',num2str(T_STEP)],'Color','r');
    text(100,140,['Images  :',fileName],'Color','r','Interpreter','none');

    %save figures
    hgsave(h_prot_edge,fullfile(dir_w, 'figures', 'prot_edge.fig'));
    print(h_prot_edge, fullfile(dir_w, 'figures', 'prot_edge.eps'),'-depsc2','-tiff');
    print(h_prot_edge, fullfile(dir_w, 'figures', 'prot_edge.tif'),'-dtiff');

    hgsave(h_prot_edge_w,fullfile(dir_w, 'figures', 'prot_edge_w.fig'));
    print(h_prot_edge_w, fullfile(dir_w, 'figures', 'prot_edge_w.eps'),'-depsc2','-tiff');
    print(h_prot_edge_w, fullfile(dir_w, 'figures', 'prot_edge_w.tif'),'-dtiff');

    if CONTR
        % plot the edge spline and the protrusion
        figure;
        hold on
        for ii=1:index
            knots=1:edge_sp_array_x(ii).number-2;
            x=fnval(edge_sp_array_x(ii),knots);
            y=fnval(edge_sp_array_y(ii),knots);
            plot(x,y);
            plot(pixel_edge{ii}(:,1),pixel_edge{ii}(:,2),'rs','MarkerSize',1);
        end
        for ii=2:index;
            quiver(protrusion{index-1}(:,1), protrusion{index-1}(:,2), protrusion{index-1}(:,3), protrusion{index-1}(:,4), 0,'r');
        end
        title('The displacement and the edge spline');
        axis([0 m_img 0  n_img]);
        axis ij;
    end
end

% Spline based version
figure(h_sup_img_edge_evolution_spline);
set(h_sup_img_edge_evolution_spline,'Visible','Off');
hold on
if PROTRUSION
    for ii=2:index
        quiver(protrusion{index-1}(:,1), protrusion{index-1}(:,2), protrusion{index-1}(:,3), protrusion{index-1}(:,4), 0,'r');
    end
end
title('Cell leading edge and protrusion vectors');
text(100,100,['Number of images  :',num2str(index)],'Color','r');
text(100,120,['Image increment   :',num2str(T_STEP)],'Color','r');
text(100,140,['Images  :',fileName],'Color','r','Interpreter','none');
set(h_sup_img_edge_evolution_spline,'Visible','On');
hgsave(h_sup_img_edge_evolution_spline,fullfile(dir_w, 'figures', 'prot_edge_spline.fig'));
print(h_sup_img_edge_evolution_spline, fullfile(dir_w, 'figures', 'prot_edge_spline.eps'),'-depsc2','-tiff');
print(h_sup_img_edge_evolution_spline, fullfile(dir_w, 'figures', 'prot_edge_spline.tif'),'-dtiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Save the overlay of speckle image and edge evolution %%
imwrite(sup_img_edge_evolution,fullfile(dir_w, 'figures', 'edge_evolution_overlay.tif'),'tif');
h_sup_img_edge_evolution = figure;
imshow(sup_img_edge_evolution);
print(h_sup_img_edge_evolution, fullfile(dir_w, 'figures', 'edge_evolution_overlay.eps'),'-depsc2','-tiff');
hgsave(h_sup_img_edge_evolution, fullfile(dir_w, 'figures', 'edge_evolution_overlay.fig'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% End plot superposition of protrusion and edge %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Calculate the assembled area per time step %%%%%%%%%%%
if MAX_IMG > 1
    index=0;
    for time=FIRST_IMG: T_STEP: FIRST_IMG + T_STEP*(MAX_IMG-1)
        index = index+1;
        %read the BW image
        fileName=char(filelist(time));
        [tmp_path,tmp_fname] = fileparts(fileName);
        mask_filename = fullfile(dir_w, 'cell_mask', ['mask_', tmp_fname, '.tif']);
        bw_img = imread(mask_filename);

        if index > 1
            %calculate current area
            [i,j,c_ta]  = find(bw_img == 1);
            cell_area(index)  = sum(sum(c_ta));

            %calculate assembly area based on BW images
            protrusion_img = 2 .* bw_img +  bw_img_old;

            %calculate created area:
            [i,j,c_a]  = find(protrusion_img == 2);
            created_area(index)   = sum(sum(c_a));
            %calculate destroid area
            [i,j,d_a] = find(protrusion_img == 1);
            destroied_area(index) = sum(sum(d_a));
        end
        %store image for next time step
        bw_img_old = bw_img;
    end

    x_array = FIRST_IMG: T_STEP: FIRST_IMG + T_STEP*(MAX_IMG-1);
    % eliminate first value
    x_array(1) = [];
    cell_area(1) = [];
    created_area(1) = [];
    destroied_area(1) = [];

    h_total_area = figure;
    plot(x_array, cell_area*PIXEL*PIXEL);
    title('Cell area');
    xlabel('Image frame #');
    ylabel('Area [nm^2]');
    box off;
    h_axes_area(1) = gca;
    h_axes_area(2) = axes('position',get(h_axes_area(1),'position'));
    xlimits = get(h_axes_area(1),'XLim');
    ylimits = get(h_axes_area(1),'YLim');
    set(h_axes_area(2),'XAxisLocation','top','YAxisLocation','right','color','none',...
        'xgrid','off','ygrid','off','box','off');
    xmin = xlimits(1) * TIME_INTERVAL;
    xmax = xlimits(2) * TIME_INTERVAL;
    ymin = ylimits(1)/PIXEL/PIXEL;
    ymax = ylimits(2)/PIXEL/PIXEL;
    axis(h_axes_area(2),[xmin xmax ymin ymax]);
    xlabel('Time step [s]');
    ylabel('Area [pixel^2]');
    hold off

    h_area = figure;
    plot(x_array, created_area*PIXEL*PIXEL);
    hold on
    plot(x_array, destroied_area*PIXEL*PIXEL,'--r');
    title('Area assembled by the cell');
    xlabel('Image frame #');
    ylabel('Area [nm^2]');
    legend('Area assembled', 'Area dissasembled');
    box off;
    h_axes_area(1) = gca;
    h_axes_area(2) = axes('position',get(h_axes_area(1),'position'));
    xlimits = get(h_axes_area(1),'XLim');
    ylimits = get(h_axes_area(1),'YLim');
    set(h_axes_area(2),'XAxisLocation','top','YAxisLocation','right','color','none',...
        'xgrid','off','ygrid','off','box','off');
    xmin = xlimits(1) * TIME_INTERVAL;
    xmax = xlimits(2) * TIME_INTERVAL;
    ymin = ylimits(1)/PIXEL/PIXEL;
    ymax = ylimits(2)/PIXEL/PIXEL;
    axis(h_axes_area(2),[xmin xmax ymin ymax]);
    %set(h_axes_area(1),'XTick',[xlimits(1):3:xlimits(2)]);
    xlabel('Time step [s]');
    ylabel('Area [pixel^2]');
    hold off

    %safe figures
    hgsave(h_total_area, fullfile(dir_w, 'figures', 'cell_area.fig'));
    print(h_total_area,  fullfile(dir_w, 'figures', 'cell_area.eps'),'-depsc2','-tiff');
    print(h_total_area,  fullfile(dir_w, 'figures', 'cell_area.tif'),'-dtiff');

    hgsave(h_area,fullfile(dir_w, 'figures', 'diff_area.fig'));
    print(h_area, fullfile(dir_w, 'figures', 'diff_area.eps'),'-depsc2','-tiff');
    print(h_area, fullfile(dir_w, 'figures', 'diff_area.tif'),'-dtiff');
    clear x_array;
end % MAX_IMG >1
%%%%%%%%%%%%%%%%%%%  End calculate the assembled area per time step %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%   Plot edge length  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_edge_length = figure;
x_array = FIRST_IMG: T_STEP: FIRST_IMG + T_STEP*(MAX_IMG-1);
plot(x_array, edge_l .* PIXEL / 1000);
title('Length of the leading edge');
xlabel('Image frame #');
ylabel('Length [\mu m]');
box off;
h_axes_edge(1) = gca;
h_axes_edge(2) = axes('position',get(h_axes_edge(1),'position'));
xlimits = get(h_axes_edge(1),'XLim');
ylimits = get(h_axes_edge(1),'YLim');
set(h_axes_edge(2),'XAxisLocation','top','YAxisLocation','right','color','none',...
    'xgrid','off','ygrid','off','box','off');
xmin = xlimits(1) * TIME_INTERVAL;
xmax = xlimits(2) * TIME_INTERVAL;
ymin = ylimits(1)/PIXEL *1000;
ymax = ylimits(2)/PIXEL *1000;
axis(h_axes_edge(2),[xmin xmax ymin ymax]);
%set(h_axes_edge(1),'XTick',[xlimits(1):3:xlimits(2)]);
xlabel('Time step [s]');
ylabel('Length [pixel]');


%safe figure
hgsave(h_edge_length,fullfile(dir_w, 'figures', 'edge_length.fig'));
print(h_edge_length, fullfile(dir_w, 'figures', 'edge_length.eps'),'-depsc2','-tiff');
print(h_edge_length, fullfile(dir_w, 'figures', 'edge_length.tif'),'-dtiff');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%   Plot control values  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(x_array, rem_pix);
h_axes_closure(1) = gca;
xlimits = get(h_axes_closure(1),'XLim');
title('Number of corrected pixels by closure operation');
xlabel('Image frame #');
ylabel('Corrected pixel');
%set(h_axes_closure(1),'XTick',[xlimits(1):3:xlimits(2)]);
clear x_array;
%%%%%%%%%%%%%%%%%   End plot control values  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Make movie  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if MOVIE
    MakeQTMovie('finish');
end
%%%%%%%%%%%%%%%%%%  End make movie  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% Averaged protrusion values  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    if PROTRUSION
        % Determine the direction of the protrusion vector
        % and store it in the assembly variable
        index=0;
        %assembly = zeros(PROTRUSION_ARRAY_SIZE, MAX_IMG);

        for time=FIRST_IMG: T_STEP: FIRST_IMG + T_STEP*(MAX_IMG-1)
            index = index +1;
            %read the BW image
            fileName=char(filelist(time));
            [tmp_path,tmp_fname] = fileparts(fileName);
            mask_filename = [dir_w  'cell_mask' filesep 'mask_' tmp_fname '.tif'];
            bw_img = imread(mask_filename);

            % non averaged protrusion vectors
            if index > 1
                for i_seg =1:nr_prot_vectors(index)
                    x = protrusion{index-1}(i_seg ,1) + protrusion{index-1}(i_seg,3);
                    y = protrusion{index-1}(i_seg ,2) + protrusion{index-1}(i_seg,4);

                    % determine if the protrusion vector is positive or negative
                    test1 = round(y) >= 1 & round(y) <= n_img;
                    test2 = round(x) >= 1 & round(x) <= m_img;
                    if test1 & test2
                        if bw_img_last(round(y), round(x)) == 1
                            assembly{index}(i_seg,index)  = -1;
                        else
                            assembly{index}(i_seg,index)  =  1;
                        end
                    end
                end %segment for
            end %index if

            %         % averaged protrusion vectors
            %         if index > 1
            %             for i_seg =1:size(av_protrusion_matrix,3)
            %                 x = av_protrusion_matrix(index,1, i_seg) + av_protrusion_matrix(index, 3, i_seg);
            %                 y = av_protrusion_matrix(index,2, i_seg) + av_protrusion_matrix(index, 4, i_seg);
            %
            %                 % determine if the protrusion vector is positive or negative
            %                 test1 = round(y) >= 1 & round(y) <= n_img;
            %                 test2 = round(x) >= 1 & round(x) <= m_img;
            %                 if test1 & test2
            %                     if bw_img_last(round(y), round(x)) == 1
            %                         av_assembly(i_seg,index)  = -1;
            %                     else
            %                         av_assembly(i_seg,index)  =  1;
            %                     end
            %                 end
            %             end %segment for
            %         end %index if

            bw_img_last =  bw_img;
        end   %time for
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %     if ~PATHS
        %         displacement_time_series = assembly .* squeeze(sqrt(disp_mech(:,1,:) .^2 + disp_mech(:,2,:) .^2));
        %     else
        %         displacement_time_series = protrusion;
        %     end

        displacement_time_series = assembly .* squeeze(sqrt(disp_mech(:,1,:) .^2 + disp_mech(:,2,:) .^2));

        %av_displacement_time_series = av_assembly' .* squeeze(sqrt(av_protrusion_matrix(:, 3, :).^2 + ...
        %                                                          av_protrusion_matrix(:, 4, :).^2));


        %find non zero entries
        displacement_non_zeros = displacement_time_series ~= 0;

        % remove first time step because there is no displacement!
        displacement_non_zeros(:,1)=[];
        displacement_time_series(:,1)=[];

        % save protrusion
        % save([dir_w  'av_abs_protrusion_matrix.mat'], 'displacement_time_series');
        %save([dir_w  'av_abs_protrusion_matrix.mat'], 'av_displacement_time_series');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%   Space average %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pos_av_protrusion = sum(displacement_time_series,1)./sum(displacement_non_zeros,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%   Time average  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        time_av_protrusion = sum(displacement_time_series,2)./sum(displacement_non_zeros,2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h_prot_velocity = figure;
        x_axis = FIRST_IMG: T_STEP: FIRST_IMG + T_STEP*(MAX_IMG-1);
        x_axis(1) = [];
        plot(x_axis,  pos_av_protrusion .* PIXEL/TIME_INTERVAL);
        box off;
        xlabel('Image frame #');
        ylabel('Protrusion velocity [nm/s]');
        h_axes_protrusion(1) = gca;
        h_axes_protrusion(2) = axes('position',get(h_axes_protrusion(1),'position'));
        xlimits = get(h_axes_protrusion(1),'XLim');
        ylimits = get(h_axes_protrusion(1),'YLim');
        set(h_axes_protrusion(2),'XAxisLocation','top','YAxisLocation','right','color','none',...
            'xgrid','off','ygrid','off','box','off');
        xmin = xlimits(1) * TIME_INTERVAL;
        xmax = xlimits(2) * TIME_INTERVAL;
        ymin = ylimits(1)/PIXEL * TIME_INTERVAL;
        ymax = ylimits(2)/PIXEL * TIME_INTERVAL;
        axis(h_axes_protrusion(2),[xmin xmax ymin ymax]);
        %set(h_axes_protrusion(1),'XTick',[xlimits(1):3:xlimits(2)]);
        xlabel('Time step [s]');
        ylabel('Displacement [pixel]');
        title('Average displacement');
        hgsave(h_prot_velocity,[dir_w  'figures' filesep 'prot_vel_pos_av.fig']);
        print(h_prot_velocity, [dir_w  'figures' filesep 'prot_vel_pos_av.eps'],'-depsc2','-tiff');
        print(h_prot_velocity, [dir_w  'figures' filesep 'prot_vel_pos_av.tif'],'-dtiff');
        hold off
        clear x_array;
        %%%%%%%%%%%%%%%%%%% End space averaged   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%  TIME AVERAGED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h_prot_velocity = figure;
        %x_axis = 1: PROT_SAMPLING : PROTRUSION_ARRAY_SIZE*PROT_SAMPLING;
        x_axis = 1: PROT_SAMPLING : size(time_av_protrusion,1)*PROT_SAMPLING;
        plot(x_axis * PIXEL /1000, time_av_protrusion .* PIXEL/TIME_INTERVAL);
        box off;
        xlabel('Position at leading edge [\mu m]');
        ylabel('Protrusion velocity [nm/s]');
        h_axes_protrusion(1) = gca;
        h_axes_protrusion(2) = axes('position',get(h_axes_protrusion(1),'position'));
        xlimits = get(h_axes_protrusion(1),'XLim');
        ylimits = get(h_axes_protrusion(1),'YLim');
        set(h_axes_protrusion(2),'XAxisLocation','top','YAxisLocation','right','color','none',...
            'xgrid','off','ygrid','off','box','off');
        xmin =  xlimits(1)/PIXEL*1000;
        xmax =  xlimits(2)/PIXEL*1000;
        ymin =  ylimits(1)/PIXEL * TIME_INTERVAL;
        ymax =  ylimits(2)/PIXEL * TIME_INTERVAL;
        axis(h_axes_protrusion(2),[xmin xmax ymin ymax]);
        %set(h_axes_protrusion(1),'XTick',[xlimits(1):1:xlimits(2)]);
        xlabel('Position at leading edge [pixel]');
        ylabel('Displacement [pixel]');
        title('Average displacement');
        hgsave(h_prot_velocity,[dir_w  'figures' filesep 'prot_vel_time_av.fig']);
        print(h_prot_velocity, [dir_w  'figures' filesep 'prot_vel_time_av.eps'],'-depsc2','-tiff');
        print(h_prot_velocity, [dir_w  'figures' filesep 'prot_vel_time_av.tif'],'-dtiff');
        hold off
        clear x_array;
        %%%%%%%%%%%%%%%%%%% End time averaged   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%% Cell activity plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if size(displacement_time_series,2) > 1
            max_prot_val = max(displacement_time_series(:));
            h_activity = figure;
            surf(displacement_time_series * PIXEL/TIME_INTERVAL,...
                'EdgeLighting','phong','EdgeColor','none','FaceColor','interp');
            title('The displacement along the leading edge ');
            xlabel('Time step');
            ylabel('Edge position [pixel]');
            zlabel('Displacement [nm/s]]');
            caxis([-max_prot_val max_prot_val]);
            colorbar;
            clear displacement_time_series
            hgsave(h_activity,[dir_w  'figures' filesep 'cell_activity.fig']);
            print(h_activity, [dir_w  'figures' filesep 'cell_activity.eps'],'-depsc2','-tiff');
            print(h_activity, [dir_w  'figures' filesep 'cell_activity.tif'],'-dtiff');
        end
        %%%%%%%%%%%%%%%%%%% End cell activity plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end %MECHANICAL & PROTRUSION
    %%%%%%%%%%%%%%%%End Averaged protrusion values  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%   Create output of the function  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~PROTRUSION
    disp_mech=0;
    img_edge=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose all
