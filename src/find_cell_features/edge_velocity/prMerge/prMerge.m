function prMerge(merg_parameters, prot_parameters)
% PRMERGE MERGES PROTRUSION, VECTORS (FLOW), SCORES (TURNOVER), IMAGES
%
%
%       This code reads the mpm.mat -containing the variables M and MPM-,
%       SCORE.mat -containing SCORE-, and the protrusion vectors. 
% 
%               SCORES has the following structure:
%               SCORES=[m,4], where
%               SCORES[m,1] is the time
%               SCORES[m,2] is the y-coordinate
%               SCORES[m,3] is the x-coordinate
%               SCORES[m,4] is the score    
%
%               M has the following structure:
%               M=[m,4,t], where
%               M[m,1,t] is the y-coordinate in t
%               M[m,2,t] is the x-coordinate in t
%               M[m,3,t] is the y-coordinate in t+1
%               M[m,4,t] is the x-coordinate in t+1    
%               
%               SCORES goes from time 2<t<T-1
%               M goes from 1<t<T-1
%
% SYNOPSIS      prMerge
%
% INPUT                :         
% 
% OUTPUT               : 
%                           
% DEPENDENCES       prMerge uses {    netAssemblyMaps
%                                     velocityMapsOrg                             
%                                       }
%
%                   prMerge is used by { 
%                                           }
%
% Matthias Machacek 12/04/03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DO_PROT         = merg_parameters.do_prot;
DO_FLOW         = merg_parameters.do_flow;
DO_SCORES       = merg_parameters.do_scores;
DO_ACTIVITY_1   = merg_parameters.do_activity_1;
DO_ACTIVITY_2   = merg_parameters.do_activity_2;

SCORES_SIGN     = merg_parameters.scores_sign;
SCORES_MATRIX   = merg_parameters.scores_matrix;
FLOW_MATRIX     = merg_parameters.flow_matrix;

% units:
% 0: sum scores in window
% 1: sum scores in window and convert it to /sec
% 2: average scores in window
% 3: average scores in window and convert it to /sec
SCORES_CONVERT  = merg_parameters.scores_convert;

% units:
% 0: pixel/frame
% 1: pixel/frame
% 2: um/min
UNITS           = merg_parameters.units;
PIXEL           = prot_parameters.pixel;
TIME_INTERVAL   = prot_parameters.time_interval;
TOTAL_TIME      = merg_parameters.total_time_steps;

SHIFT_SCORES_HALF_TIME_STEP = merg_parameters.shift_scores_half_time_step; 
% if REPLACE_NAN = 1 then all missing data points are replaced 
% with zeros
REPLACE_NAN     = merg_parameters.replace_nan; 
PROJ_FLOW       = 1;

if DO_PROT
    PROT_SAMPLING   = prot_parameters.prot_sampling;
    PARENTH_R       = prot_parameters.parenth_r;
    PARENTH_L       = prot_parameters.parenth_l;
else
    PROT_SAMPLING   = prot_parameters.prot_sampling;
    PARENTH_R       = 0;
    PARENTH_L       = 0;  
end

PROJECT_DIR     = merg_parameters.project_dir;
PROT_DIR        = merg_parameters.prot_dir;
FLOW_FILE       = merg_parameters.flow_dir;
SCORES_DIR      = merg_parameters.scores_dir;

ACTIVITY_DIR_1  = merg_parameters.activity_dir_1;
ACTIVITY_DIR_2  = merg_parameters.activity_dir_2;

SEG_NR          = merg_parameters.seg_nr;
SEG_LENGTH      = merg_parameters.seg_length;
SEG_DEPTH       = merg_parameters.seg_depth;

WINDOWS_TYPE    = merg_parameters.windows_type;

% disable this parameter!
FIRST_TIME      = merg_parameters.first_time;
TOTAL_TIME_STEPS= merg_parameters.total_time_steps;
START_SEG       = merg_parameters.start_seg;
START_SEG = START_SEG+1;
END_SEG         = merg_parameters.end_seg;

SEG_SHIFT       = merg_parameters.seg_shift;
PLOT_VECTORS    = merg_parameters.plot_vectors;
INTERPOLATE     = merg_parameters.interpolation;
D0              = merg_parameters.d0;

%Get the time window and alignment parameters.
protTimeStep     = merg_parameters.protTimeStep;
protTimeWinL     = merg_parameters.protTimeWinL;
flowTimeStep     = merg_parameters.flowTimeStep;
flowTimeWinL     = merg_parameters.flowTimeWinL;
scoreTimeStep    = merg_parameters.scoreTimeStep;
scoreTimeWinL    = merg_parameters.scoreTimeWinL;
activityTimeStep = merg_parameters.activityTimeStep;
activityTimeWinL = merg_parameters.activityTimeWinL;

% WINDOWS_TYPE == 1 dynamics windows 
% WINDOWS_TYPE == 2 static window
% WINDOWS_TYPE == 3 static window band

if WINDOWS_TYPE == 2  
    SEG_NR          = 1;
    SEG_LENGTH      = 0;
    SEG_DEPTH       = 0;
end

warning off all;

% fill the activity structure
if DO_ACTIVITY_1
    DO_ACTIVITY(1) = 1;
else
    DO_ACTIVITY(1) = 0;
end
if DO_ACTIVITY_2
    DO_ACTIVITY(2) = 1;
else
    DO_ACTIVITY(2) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Parameters for the smoothing of the activity images %%%%%%%%%%%%%%%%%%
IMG_GAUSS_W     = 20; 
IMG_GAUSS_SIG   = 6;  
IMAGE_STRETCH   = 15;

% Number of time steps used for the averaging of the scores
DELTA_T_AV = 1;

% Determine what units to use
if UNITS == 0
    % use pixel and frame
    PIXEL           = 1;
    TIME_INTERVAL   = 1;
elseif UNITS == 1
    % keep nm and s
    PIXEL           = PIXEL;
    TIME_INTERVAL   = TIME_INTERVAL;
elseif UNITS == 2
    % convert nm to um and sec to min
    PIXEL           =  PIXEL / 1000;
    TIME_INTERVAL   =  TIME_INTERVAL / 60;    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Create directories    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the directury for merged files   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
merg_dir = [PROJECT_DIR filesep 'merg'];
exist_dir = exist(merg_dir,'dir');
if exist_dir == 0
    %it does not exist, so try to create it
    [s, mess, messid] = mkdir(merg_dir);
    if s==0
        disp('Failed to create the merg directory');
        return
    end
end
merg_dir = [merg_dir filesep];


% Create directory for overlay images     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PLOT_VECTORS
    exist_dir = exist([merg_dir 'overlay']);
    if exist_dir == 0
        %it does not exist, so try to create it
        [s, mess, messid] = mkdir([merg_dir 'overlay']);
        if s==0
            disp('Failed to create overlay directory');
            return
        end
    end
end

% Create directory for merged data (.mat files)    %%%%%%%%%%%%%%%%%%%%%%%%
exist_dir = exist([merg_dir 'data']);
if exist_dir == 0
    %it does not exist, so try to create it
    [s, mess, messid] = mkdir([merg_dir 'data']);
    if s==0
        disp('Failed to create data directory');
        return
    end
end

% Create directory for averaged data (.mat files)   %%%%%%%%%%%%%%%%%%%%%%%
exist_dir = exist([merg_dir 'av_data']);
if exist_dir == 0
    %it does not exist, so try to create it
    [s, mess, messid] = mkdir([merg_dir 'av_data']);
    if s==0
        disp('Failed to create data directory');
        return
    end
end
av_data_dir = [merg_dir 'av_data' filesep];


% Create directory for figures        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist_dir == 0
    %it does not exist, so try to create it
    [s, mess, messid] = mkdir([merg_dir 'figures']);
    if s==0
        disp('Failed to create figure directory');
        return
    end
end
figures_dir = [merg_dir 'figures' filesep];

%%%%%%%%%%%   End create directories    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%    Get file names lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mask file names
filelist_temp = dir([PROJECT_DIR filesep PROT_DIR filesep 'cell_mask']);
maskImgIndex = [];
iEntry = 1;
for( i = 1:length(filelist_temp))
   if(~filelist_temp(i).isdir)
      filelist_mask(iEntry) = {[PROJECT_DIR filesep PROT_DIR filesep 'cell_mask' filesep filelist_temp(i).name]};
      [path,body,no,ext]=getFilenameBody(filelist_mask{iEntry});
      maskImgIndex = [maskImgIndex str2num(no)];
      iEntry = iEntry + 1;
   end
end
clear filelist_temp


% If we generate overlay images we have to determine which images to use
if PLOT_VECTORS
% get image file names

[IMG_DIR, name, ext, versn] = fileparts(prot_parameters.file);   
%    %[IMG_DIR, name, ext, versn] = fileparts('M:\GTPases\cdc42\Masked_data\exp1_movie1-1\activity_MEF_ratio1_123104_lp_corrected\crop_ratio_corr_default01.tif');
%    %[IMG_DIR, name, ext, versn] = fileparts('S:\scripps\analysis\machacek\FocalAdhesions\80505\set_2\pax2_cut1\crop_default001.tif');
%    
%    [name, path] = uigetfile('*','Select first image for segment overlay');
%    if path == 0
%        return;
%    end 
%    
    filelist_images_b = dir([IMG_DIR filesep '*.tif']);
%    if length(filelist_images_b) == 0
%        filelist_images_b = dir([path '*.bmp']);
%    end
%

    iEntry = 1;
    for i=1:length(filelist_images_b)
        if(~filelist_images_b(i).isdir)
            filelist_images(iEntry,:) = filelist_images_b(i).name;
            iEntry = iEntry + 1;
        end
    end
    [path,body,no,ext]=getFilenameBody(filelist_images(1,:));
    firstImgIndex = str2num(no);
end

% file name of the the pixel edge
file_pixel_edge = [PROJECT_DIR  filesep PROT_DIR filesep 'pixel_edge.mat'];

% file name of the normal vectors
file_normal=[PROJECT_DIR filesep PROT_DIR filesep 'normal_matrix.mat'];

% get the filelist filenames containing the protrusion vectors    
if DO_PROT
  	file_protrusion=[PROJECT_DIR filesep PROT_DIR filesep 'protrusion.mat'];
end

% file name containing the scores maps
if DO_SCORES
    if SCORES_MATRIX == 1
        file_scores=[SCORES_DIR];
    else
      	filelist_scores_b = dir([SCORES_DIR filesep '*']);
        iEntry = 1;
     	for i=1:length(filelist_scores_b)
          	if(~filelist_scores_b(i).isdir) 
                filelist_scores(iEntry,:) = filelist_scores_b(i).name;
                iEntry = iEntry + 1;
            end
        end     
    end
end

% get the filelist containing the velocity vectors
if DO_FLOW
 	if FLOW_MATRIX == 1
        file_vectors= FLOW_FILE;
    else
        [FLOW_PATH, name, ext, versn] = fileparts(FLOW_FILE);
        filelist_flow_b = dir([FLOW_PATH filesep '*']);
        iEntry = 1;
    	for i=1:length(filelist_flow_b)
            if(~filelist_flow_b(i).isdir)
                filelist_flow(iEntry,:) = filelist_flow_b(i).name;
                iEntry = iEntry + 1;
            end
        end         
    end   
end


% file name containing the activity images
if DO_ACTIVITY_1
    filelist_activity_b = dir([ACTIVITY_DIR_1 filesep '*.tif']);
    % if there are none check for bmp format images
    if length(filelist_activity_b) == 0
        filelist_activity_b = dir([ACTIVITY_DIR_1 filesep '*.bmp']);
    end
    for i=1:length(filelist_activity_b)
        filelist_activity_1(i,:) = filelist_activity_b(i).name;
    end
end
if DO_ACTIVITY_2
    % file name containing the activity images
    filelist_activity_b = dir([ACTIVITY_DIR_2 filesep '*.tif']);
    if length(filelist_activity_b) == 0
        filelist_activity_b = dir([ACTIVITY_DIR_2 filesep '*.bmp']);
    end
    for i=1:length(filelist_activity_b)
        filelist_activity_2(i,:) = filelist_activity_b(i).name;
    end
end
%%%%%%%%%%%%%% End get file name list  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Write all the parameters to a file  %%%%%%%%%%%%%%%%%%%%%%%%%%%
save_pralpha_parameters([merg_dir 'merg_parameters.dat'], merg_parameters);
save([merg_dir  'merg_parameters'], 'merg_parameters');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% load data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = load(file_pixel_edge);
pixel_edge = s.pixel_edge;

% load data with edge splines
s = load([PROJECT_DIR filesep  PROT_DIR filesep 'edge_spline.mat']);
edge_sp_array_x = s.edge_sp_array_x;
edge_sp_array_y = s.edge_sp_array_y;

% this loads cell normal_matrix 
s = load(file_normal);
normal_matrix = s.normal_matrix;

if DO_PROT
    % this loads cell 'protrusion'
    load(file_protrusion);
end

if DO_SCORES && SCORES_MATRIX
    % read the poly/depoly score map
    scores_var = load(file_scores);
    a = fieldnames(scores_var);
    scores_name = char(a);
    scores_matrix = scores_var.(scores_name);
    clear scores_var;
end

if DO_FLOW && FLOW_MATRIX
    % read the velocity vectors
    flow_var = load(file_vectors);
    a = fieldnames(flow_var);
    if length(a) == 2
        % there is a M and MPM variable in flow_var
        flow_name = char(a{1});
        flow = flow_var.(flow_name);       
    else
        flow_name = char(a);
        flow = flow_var.(flow_name);
    end
    clear flow_var;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if WINDOWS_TYPE == 2

    % get the static window from the user
    ans = questdlg('Use existing region?','Just do it','Use existing','Define new','I dont know');
    
    if strcmp(ans, 'Define new');
        helpdlg('Select image to indicate window');
        [name, path] = uigetfile('*','Select image to indicate window');
        if path == 0
            return;
        end
        img_window = imread([path name]);
        img_window = imadjust(img_window, stretchlim(img_window),[]);
        helpdlg('Select region of interested to collect data');
        figure
        [static_window,xb,yb]  = roipoly(img_window);
        
        % save it 
        save([merg_dir sampling_window.mat'], 'static_window','xb','yb');
        
        % get the corresponding protrusion region from the user
        img_window = imadjust(img_window, stretchlim(img_window),[]);
        helpdlg('Select connected edge region of interested for protrusion');
        figure
        static_window_prot = roipoly(img_window);
        save([merg_dir prot_sampling_window.mat'],'static_window_prot');
    else
        % let user specify
        [fn_tmp,pn_tmp] = uigetfile('FilterSpec','Select variable region');
        % static_window 
        load([pn_tmp, fn_tmp]);
        
        [fn_tmp,pn_tmp] = uigetfile('FilterSpec','Select protrusion region');
        % static_window_prot
        load([pn_tmp, fn_tmp]);
    end
elseif  WINDOWS_TYPE == 3
    
    
    
    
end

%the first time step refers to number of the files in the folder% 
%and NOT to the number in the file name!!!!!!!!!!!!   %%%%%%%%%%%
time_index = 0;
strg = sprintf('%%.%dd',3);

for time = FIRST_TIME : 1: TOTAL_TIME_STEPS
%for time = 1 : 1: TOTAL_TIME_STEPS
    time_index = time_index+1;
    
    backSpc = ['\b\b\b'];
    fprintf(1,[strg],time_index);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the pixel edge   %%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_p_edge  = pixel_edge{time}(:,1);
    y_p_edge  = pixel_edge{time}(:,2);
    n_pix     = size(x_p_edge,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Determine the length of a segment (pixel) %%%%%%%%%%%%%
    dl = 0;
    i_l = PARENTH_L * PROT_SAMPLING;
    i_r = PARENTH_R * PROT_SAMPLING;
    for i=i_l+1 : n_pix -1-i_r
        dl = dl+sqrt((x_p_edge(i+1) - x_p_edge(i))^2+(y_p_edge(i+1) - y_p_edge(i))^2);
    end
    segment_length(time_index) = dl/SEG_NR;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fileName_mask=char(filelist_mask(time));
    mask_img=imread(fileName_mask);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%% Get image size from the image mask  %%%%%%%%%%%%%%%%%%
    if time_index == 1
        imgSize = size(mask_img);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the activities  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DO_ACTIVITY_1
        fileName_activity_1=[ACTIVITY_DIR_1 filesep char(filelist_activity_1(time,:))];
        activity_map{1} = imread(fileName_activity_1);
    end
    if DO_ACTIVITY_2
        fileName_activity_2=[ACTIVITY_DIR_2 filesep char(filelist_activity_2(time,:))];
        activity_map{2} = imread(fileName_activity_2);
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  average the unit normals in segments %%%%%%%%%%%%%%%%%%%%%%
    p_n=1:3:n_pix; % the number 3 comes from the imEdgeTracker.m
    if WINDOWS_TYPE == 1
        % average the unit normals in segments
       
        [x_normal_out_av, y_normal_out_av, x_av_pos_normal, y_av_pos_normal, m_pos]=...
            prGetAvEdge(edge_sp_array_x(time), edge_sp_array_y(time), edge_sp_array_x(time).knots(end),...
            p_n, normal_matrix{time}(:,3), normal_matrix{time}(:,4), 'nr_sect', SEG_NR);

        % re-normalize the averaged normals
        l_n_av = sqrt(x_normal_out_av.^2 + y_normal_out_av.^2);
        x_av_normal = x_normal_out_av ./ l_n_av;
        y_av_normal = y_normal_out_av ./ l_n_av;
    elseif WINDOWS_TYPE == 2
        x = int16(fnval(edge_sp_array_x(time), 1));
        y = int16(fnval(edge_sp_array_y(time), 1));
        i=1;
        while p_n(i) < n_pix && ~static_window_prot(y,x) 
            x = int16(fnval(edge_sp_array_x(time), p_n(i)));
            y = int16(fnval(edge_sp_array_y(time), p_n(i)));
            i=i+1;
        end
        % now accumulate
        x_normal_out_av = 0;
        y_normal_out_av = 0;
        x_av_pos_normal = 0;
        y_av_pos_normal = 0;
        num_el = 0;
        while p_n(i) < n_pix && static_window_prot(int16(y),int16(x)) 
            x = fnval(edge_sp_array_x(time), p_n(i));
            y = fnval(edge_sp_array_y(time), p_n(i));
            
            x_av_pos_normal = x_av_pos_normal + x;
            y_av_pos_normal = y_av_pos_normal + y;
             
            x_normal_out_av = x_normal_out_av+normal_matrix{time}(i,3);
            y_normal_out_av = y_normal_out_av+normal_matrix{time}(i,4);
             
            num_el = num_el + 1;
            i=i+1;          
        end        
        x_normal_out_av = x_normal_out_av / num_el;
        y_normal_out_av = y_normal_out_av / num_el;
        x_av_pos_normal = x_av_pos_normal / num_el;
        y_av_pos_normal = y_av_pos_normal / num_el;
        
        % re-normalize the averaged normals
        l_n_av = sqrt(x_normal_out_av.^2 + y_normal_out_av.^2);
        x_av_normal = x_normal_out_av ./ l_n_av;
        y_av_normal = y_normal_out_av ./ l_n_av;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  average the protrusion in segments %%%%%%%%%%%%%%%%%%%%%%%%
    if DO_PROT
        % Average the protrusion vectors in the segments
        i_nn=1:PROT_SAMPLING:n_pix;
        l_i_nn=length(i_nn);
        i_nn(l_i_nn-PARENTH_R:l_i_nn)=[];
        i_nn(1:PARENTH_L)=[];
        clear l_i_nn;
        if WINDOWS_TYPE == 1
            [x_av_prot, y_av_prot, x_av_pos_prot, y_av_pos_prot, spline_p_av_prot]=...
                prGetAvEdge(edge_sp_array_x(time), edge_sp_array_y(time), n_pix,...
                i_nn, protrusion{time}(:,3), protrusion{time}(:,4), 'nr_sect', SEG_NR);

        elseif WINDOWS_TYPE == 2
            x = int16(fnval(edge_sp_array_x(time), 1));
            y = int16(fnval(edge_sp_array_y(time), 1));
            i=1;
            while i_nn(i) < n_pix && ~static_window_prot(y,x)
                x = int16(fnval(edge_sp_array_x(time), i_nn(i)));
                y = int16(fnval(edge_sp_array_y(time), i_nn(i)));
                i=i+1;
            end
            % now accumulate
            x_av_prot = 0;
            y_av_prot = 0;
            x_av_pos_prot = 0;
            y_av_pos_prot = 0;
            num_el = 0;
            while i_nn(i) < n_pix && static_window_prot(int16(y),int16(x))
                x = fnval(edge_sp_array_x(time), i_nn(i));
                y = fnval(edge_sp_array_y(time), i_nn(i));

                x_av_pos_prot = x_av_pos_prot + x;
                y_av_pos_prot = y_av_pos_prot + y;

                x_av_prot = x_av_prot+protrusion{time}(i,3);
                y_av_prot = y_av_prot+protrusion{time}(i,4);
                num_el = num_el +1;
                i=i+1;
            end
            x_av_prot = x_av_prot / num_el;
            y_av_prot = y_av_prot / num_el;
            x_av_pos_prot = x_av_pos_prot / num_el;
            y_av_pos_prot = y_av_pos_prot / num_el;
        end
    end
    %n_el_av_prot = SEG_NR + 1;
    %n_el_av_norm = SEG_NR + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% get filtered score maps %%%%%%%%%%%%%%%%%%%%%%%%%
    if DO_SCORES
     	if SCORES_MATRIX
            %[score_map, score_map_sign, num_t_scores(time_index)] =...
            %      netAssemblyMaps(scores, time, DELTA_T_AV, imgSize, contr);  
            score_data = scores_matrix(find(scores_matrix(:,1) == time),:);
        else
            fileName_scores=[SCORES_DIR filesep char(filelist_scores(time,:))]; 
            score_tmp = load(fileName_scores);
            a=char(fieldnames(score_tmp));
            score_data = score_tmp.(a);
            clear score_tmp;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Read the velocity map %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DO_FLOW
        if FLOW_MATRIX
            velocity_field = prExtractVelocities(imgSize, flow, time, mask_img, D0);
        else
            fileName_flow=[FLOW_PATH filesep char(filelist_flow(time,:))]; 
            flow = load(fileName_flow);
            
            a=char(fieldnames(flow));
            velocity_field = flow.(a);
            
            %Added by Lin Ji (Jan 29, 2007). Remove 'NaN'.
            nanInd = find(isnan(velocity_field(:,1)) | isnan(velocity_field(:,2)) | ...
               isnan(velocity_field(:,3)) | isnan(velocity_field(:,4)));
            velocity_field(nanInd,:) = [];
            
            velocity_field(:,3) = velocity_field(:,3) - velocity_field(:,1);
            velocity_field(:,4) = velocity_field(:,4) - velocity_field(:,2);
            % remove Nan's
            nan_ind = find(isnan(velocity_field(:,3)));
            velocity_field(nan_ind,:) = [];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
   
    
    seg_index=0;
    for seg = START_SEG: SEG_NR - END_SEG
        seg_index=seg_index+1;
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Generate segments  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if WINDOWS_TYPE == 1
            % use dynamics windows
            [xb, yb] = prCreateSegments(edge_sp_array_x(time), edge_sp_array_y(time),...
                x_av_normal, y_av_normal, imgSize(2), imgSize(1), SEG_NR, SEG_DEPTH, seg);
            %%%%%%%%%%%%%%%%% Create the same segment placed some way %%%%%%%%
            %%%%%%%%%%%%%%%%% Set back along the normal to the edge %%%%%%%%%%
            if strcmp(merg_parameters.segShiftDir,'normal')
              xb_shift  = xb - SEG_SHIFT * x_av_normal(seg);
              yb_shift  = yb - SEG_SHIFT * y_av_normal(seg);
            else
            %%%%%%%%% Set back parallel to average normal of segment normals %%%%%%%%%%
              xb_shift  = xb - SEG_SHIFT * mean(x_av_normal);
              yb_shift  = yb - SEG_SHIFT * mean(y_av_normal);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Special Lin feature: shift leading edge part of segments forward
            if 0
                forward_shift = 10;   % in pixel
                xb  = xb + forward_shift * x_av_normal(seg);
                yb  = yb + forward_shift * y_av_normal(seg);
            end
        elseif WINDOWS_TYPE == 2
            % use single static window
            %[xb, yb] = static_window;


        end
        %%%%%%%%%%%% End Generate segments  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% PROTRUSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        if DO_PROT
            av_prot         = sqrt(x_av_prot(seg)^2 + y_av_prot(seg)^2);

            % Protrusion into the normal direction
            scal            = x_av_normal(seg) .* x_av_prot(seg) + y_av_normal(seg) .* y_av_prot(seg);
            x_av_prot_proj  = scal * x_av_normal(seg);
            y_av_prot_proj  = scal * y_av_normal(seg);
            av_prot_proj    = sign(scal)*sqrt(x_av_prot_proj^2 + y_av_prot_proj^2);

            % get the angle between the normal direction and the retrograde flow
            % consider just the absolute difference without the direction !!
            protrusion_angle(seg_index,time_index)   = 180/pi * real(acos(abs(scal) / av_prot));
            % rescale
            protrusion_seg(seg_index,time_index)         = av_prot     * PIXEL / TIME_INTERVAL;
            protrusion_normal(seg_index,time_index)      = av_prot_proj* PIXEL / TIME_INTERVAL;
            %protrusion_normal(seg_index,time_index)=sign(scal)*protrusion(
            %seg_index,time_index);
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%  FLOW  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DO_FLOW
            % Extract from the velocity map
            velocity_field_crop = velocity_field;
            index=inpolygon(velocity_field(:,1),velocity_field(:,2),yb,xb);
            velocity_field_crop(find(~index),:)=[];

            % Average retrograde flow
            if size(velocity_field_crop,1) > 0
                % Number of elements in the segment
                retrograde_flow_num(seg_index,time_index) = size(velocity_field_crop,1);

                % Average the velocity field in segment and get the variance
                av_velocity_field_crop          = sum(velocity_field_crop,1) ./ retrograde_flow_num(seg_index,time_index);
                var_seg_velocities              = var(velocity_field_crop(:,4));

                
                av_vector_field(seg_index,time_index,1) = av_velocity_field_crop(3);
                av_vector_field(seg_index,time_index,2) = av_velocity_field_crop(4);
                
                
                av_retrograde_flow_crop         = sqrt(av_velocity_field_crop(4)^2 + av_velocity_field_crop(3)^2);
                % Project the average retrograde flow onto the normal direction
                scal                            = x_av_normal(seg) .* av_velocity_field_crop(:,4) + ...
                                                  y_av_normal(seg) .* av_velocity_field_crop(:,3);
                if PROJ_FLOW                           
                    x_av_velocity_field_crop_proj   = scal * x_av_normal(seg);
                    y_av_velocity_field_crop_proj   = scal * y_av_normal(seg);
                    av_retrograde_flow_crop_proj    = sign(scal)* sqrt(x_av_velocity_field_crop_proj^2 + ...
                                                                       y_av_velocity_field_crop_proj^2);
                else
                    av_retrograde_flow_crop_proj    = av_retrograde_flow_crop; 
                end
                % Get the angle between the normal direction and the retrograde flow
                % here, because the normal can be slightly longer then one (eg.0005) take the real part
                retrograde_flow_angle(seg_index,time_index)     = 180/pi * real(acos(abs(scal) / av_retrograde_flow_crop));
                % Rescale!
                retrograde_flow(seg_index,time_index)           = av_retrograde_flow_crop      * PIXEL / TIME_INTERVAL;
                retrograde_flow_normal(seg_index,time_index)    = av_retrograde_flow_crop_proj * PIXEL / TIME_INTERVAL;
            else
                av_velocity_field_crop                          = [0 0 0 0];
                retrograde_flow_num(seg_index,time_index)       = 0;
                retrograde_flow_angle(seg_index,time_index)     = NaN;
                retrograde_flow(seg_index,time_index)           = NaN;
                retrograde_flow_normal(seg_index,time_index)    = NaN;
                av_vector_field(seg_index,time_index,1)         = NaN;
                av_vector_field(seg_index,time_index,2)         = NaN;
            end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%% SHIFTED RETROGRADE FLOW %%%%%%%%%%%%%%%%%%%%%%%%%
            %% Extract from the velocity map
            if WINDOWS_TYPE == 1
            shift_velocity_field_crop  = velocity_field;
            index  = inpolygon(velocity_field(:,1), velocity_field(:,2), yb_shift,  xb_shift);
            shift_velocity_field_crop(find(~index),:)  =[];


            % Average shifted retrograde flow
            if size(shift_velocity_field_crop,1) > 0
                % Number of elements in the segment
                shift_retrograde_flow_num(seg_index,time_index) = size(shift_velocity_field_crop,1);

                % Average the velocity field in segment and calc the variance
                av_shift_velocity_field_crop = sum(shift_velocity_field_crop,1) ./ shift_retrograde_flow_num(seg_index,time_index);
                %var_seg_shift_velocities(seg_index)=var(shift_velocity_field_crop(:,4));

                av_shift_vector_field(seg_index,time_index,1) = av_shift_velocity_field_crop(3);
                av_shift_vector_field(seg_index,time_index,2) = av_shift_velocity_field_crop(4);
                
                
                av_shift_retrograde_flow_crop           = sqrt(av_shift_velocity_field_crop(4)^2 + av_shift_velocity_field_crop(3)^2) ;
                % average shifted retrograde flow in the normal direction
                scal                                    = x_av_normal(seg) .* av_shift_velocity_field_crop(:,4) + y_av_normal(seg) .* av_shift_velocity_field_crop(:,3);
                x_av_shift_velocity_field_crop_proj     = scal * x_av_normal(seg);
                y_av_shift_velocity_field_crop_proj     = scal * y_av_normal(seg);
                av_shift_retrograde_flow_crop_proj      = sign(scal)* sqrt(x_av_shift_velocity_field_crop_proj^2+y_av_shift_velocity_field_crop_proj^2);
                % Rescale!
                shift_retrograde_flow(seg_index,time_index)         = av_shift_retrograde_flow_crop      * PIXEL / TIME_INTERVAL;
                shift_retrograde_flow_normal(seg_index,time_index)  = av_shift_retrograde_flow_crop_proj * PIXEL / TIME_INTERVAL;
            else
                av_shift_velocity_field_crop = [0 0 0 0];
                shift_retrograde_flow_num(seg_index,time_index)     = 0;
                shift_retrograde_flow_normal(seg_index,time_index)  = NaN;
                av_shift_vector_field(seg_index,time_index,1)       = NaN;
                av_shift_vector_field(seg_index,time_index,2)       = NaN;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compare the direction of network flow and shifted network
            % flow
            av_network_velocity_direction_cor(seg_index,time_index) = ...
                av_velocity_field_crop(3) * av_shift_velocity_field_crop(3)+...
                av_velocity_field_crop(4) * av_shift_velocity_field_crop(4);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end % WINDOWS_TYPE == 1
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% SCORES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DO_SCORES
            % Extract the values
            clear in;
            in = inpolygon(score_data(:,3),score_data(:,2),xb,yb);
            score_data_extr = [score_data(in,3), score_data(in,2), score_data(in,4)];
            score_data_extr_plus = [score_data(in,3), score_data(in,2), (score_data(in,4)>0).*score_data(in,4)];
            score_data_extr_minus = [score_data(in,3), score_data(in,2), (score_data(in,4)<0).*score_data(in,4)];
            
            % Get number of scores in segment
            score_num(seg_index,time_index) = size(score_data_extr,1);
            score_num_plus(seg_index,time_index) = sum(score_data_extr_plus(1:size(score_data_extr_plus),3)>0);
            score_num_minus(seg_index,time_index) = sum(score_data_extr_minus(1:size(score_data_extr_minus),3)<0);
%             score_num_plus(seg_index,time_index) = size(score_data_extr_plus,1);
%             score_num_minus(seg_index,time_index) = size(score_data_extr_minus,1);
            
            if score_num(seg_index,time_index) > 0
                av_score                        = sum(score_data_extr(:,3));
                % Rescale!% units:
                % 0: sum scores in window
                % 1: sum scores in window and convert it to /sec
                % 2: average scores in window
                % 3: average scores in window and convert it to /sec
                if SCORES_CONVERT == 0
                    score(seg_index,time_index)     = av_score;
                elseif SCORES_CONVERT == 1
                    score(seg_index,time_index)     = av_score / TIME_INTERVAL;
                elseif SCORES_CONVERT == 2
                    score(seg_index,time_index)     = av_score / score_num(seg_index,time_index);                    
                elseif SCORES_CONVERT == 3
                    score(seg_index,time_index)     = av_score /...
                        score_num(seg_index,time_index) / TIME_INTERVAL;         
                end
            else
                score(seg_index,time_index)     = NaN;
                score_num(seg_index,time_index) = 0;
            end
            if score_num_plus(seg_index,time_index) > 0
                av_score_plus                        = sum(score_data_extr_plus(:,3));
                % Rescale!% units:
                % 0: sum scores_plus in window
                % 1: sum scores_plus in window and convert it to /sec
                % 2: average scores_plus in window
                % 3: average scores_plus in window and convert it to /sec
                if SCORES_CONVERT == 0
                    score_plus(seg_index,time_index)     = av_score_plus;
                elseif SCORES_CONVERT == 1
                    score_plus(seg_index,time_index)     = av_score_plus / TIME_INTERVAL;
                elseif SCORES_CONVERT == 2
                    score_plus(seg_index,time_index)     = av_score_plus / score_num_plus(seg_index,time_index);                    
                elseif SCORES_CONVERT == 3
                    score_plus(seg_index,time_index)     = av_score_plus /...
                        score_num_plus(seg_index,time_index) / TIME_INTERVAL;         
                end
            else
                score_plus(seg_index,time_index)     = NaN;
                score_num_plus(seg_index,time_index) = 0;
            end
            if score_num_minus(seg_index,time_index) > 0
                av_score_minus                        = sum(score_data_extr_minus(:,3));
                % Rescale!% units:
                % 0: sum scores_plus in window
                % 1: sum scores_plus in window and convert it to /sec
                % 2: average scores_plus in window
                % 3: average scores_plus in window and convert it to /sec
                if SCORES_CONVERT == 0
                    score_minus(seg_index,time_index)     = av_score_minus;
                elseif SCORES_CONVERT == 1
                    score_minus(seg_index,time_index)     = av_score_minus / TIME_INTERVAL;
                elseif SCORES_CONVERT == 2
                    score_minus(seg_index,time_index)     = av_score_minus / score_num_minus(seg_index,time_index);                    
                elseif SCORES_CONVERT == 3
                    score_minus(seg_index,time_index)     = av_score_minus /...
                        score_num_minus(seg_index,time_index) / TIME_INTERVAL;         
                end
            else
                score_minus(seg_index,time_index)     = NaN;
                score_num_minus(seg_index,time_index) = 0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% SHIFTED SCORES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if WINDOWS_TYPE == 1
        if DO_SCORES       
             % Extract the values
             clear in_shift;
             in_shift = inpolygon(score_data(:,3),score_data(:,2),xb_shift,yb_shift);
             score_data_extr = [score_data(in_shift,3), score_data(in_shift,2), score_data(in_shift,4)];
             score_data_extr_plus = [score_data(in_shift,3), score_data(in_shift,2), (score_data(in_shift,4)>0).*score_data(in_shift,4)];
             score_data_extr_minus = [score_data(in_shift,3), score_data(in_shift,2), (score_data(in_shift,4)<0).*score_data(in_shift,4)];
            % Get number of scores in shifted segment
            shift_score_num(seg_index,time_index) = size(score_data_extr,1);
            shift_score_num_plus(seg_index,time_index) = sum(score_data_extr_plus(1:size(score_data_extr_plus,1),3)>0);
            shift_score_num_minus(seg_index,time_index) = sum(score_data_extr_minus(1:size(score_data_extr_minus,1),3)<0);
%             shift_score_num_plus(seg_index,time_index) = size(score_data_extr_plus,1);
%             shift_score_num_minus(seg_index,time_index) = size(score_data_extr_minus,1);
            
            if shift_score_num(seg_index,time_index) > 0
                av_score                        = sum(score_data_extr(:,3));
                % Rescale! This is score per second!
                if SCORES_CONVERT == 0
                    shift_score(seg_index,time_index)     = av_score;
                elseif SCORES_CONVERT == 1
                    shift_score(seg_index,time_index)     = av_score / TIME_INTERVAL;
                elseif SCORES_CONVERT == 2
                    shift_score(seg_index,time_index)     = av_score / shift_score_num(seg_index,time_index);                    
                elseif SCORES_CONVERT == 3
                    shift_score(seg_index,time_index)     = av_score /...
                        shift_score_num(seg_index,time_index) / TIME_INTERVAL;         
                end
            else
                shift_score(seg_index,time_index)     = NaN;
                shift_score_num(seg_index,time_index) = 0;
            end 
            
           if shift_score_num_plus(seg_index,time_index) > 0
                av_score_plus                        = sum(score_data_extr_plus(:,3));
                % Rescale! This is score per second!
                if SCORES_CONVERT == 0
                    shift_score_plus(seg_index,time_index)     = av_score_plus;
                elseif SCORES_CONVERT == 1
                    shift_score_plus(seg_index,time_index)     = av_score_plus / TIME_INTERVAL;
                elseif SCORES_CONVERT == 2
                    shift_score_plus(seg_index,time_index)     = av_score_plus / shift_score_num_plus(seg_index,time_index);                    
                elseif SCORES_CONVERT == 3
                    shift_score_plus(seg_index,time_index)     = av_score_plus /...
                        shift_score_num_plus(seg_index,time_index) / TIME_INTERVAL;         
                end
            else
                shift_score_plus(seg_index,time_index)     = NaN;
                shift_score_num_plus(seg_index,time_index) = 0;
           end 
           
            if shift_score_num_minus(seg_index,time_index) > 0
                av_score_minus                        = sum(score_data_extr_minus(:,3));
                % Rescale! This is score per second!
                if SCORES_CONVERT == 0
                    shift_score_minus(seg_index,time_index)     = av_score_minus;
                elseif SCORES_CONVERT == 1
                    shift_score_minus(seg_index,time_index)     = av_score_minus / TIME_INTERVAL;
                elseif SCORES_CONVERT == 2
                    shift_score_minus(seg_index,time_index)     = av_score_minus / shift_score_num_minus(seg_index,time_index);                    
                elseif SCORES_CONVERT == 3
                    shift_score_minus(seg_index,time_index)     = av_score_minus /...
                        shift_score_num_minus(seg_index,time_index) / TIME_INTERVAL;         
                end
            else
                shift_score_minus(seg_index,time_index)     = NaN;
                shift_score_num_minus(seg_index,time_index) = 0;
           end 
        end
        end % WINDOWS_TYPE == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% ACTIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for activity_i = 1:2
            if DO_ACTIVITY(activity_i)
                % Create mask
                img=zeros(imgSize);
                mask_image=roipoly(img,xb,yb);

                % Extract from the activity map
                activity_map_crop = double(mask_image) .* double(activity_map{activity_i});

                clear activity_map_crop_extr;
                % Find the scores in the extracted segment of the score map
                [activity_map_crop_extr(:,2), activity_map_crop_extr(:,1), activity_map_crop_extr(:,3)] = find(activity_map_crop);

                % Get number of scores in segment
                activity_num(activity_i, seg_index,time_index) = size(activity_map_crop_extr,1);

                % Average activity map
                if activity_num(activity_i, seg_index,time_index) > 0
                    av_activity                        = sum(activity_map_crop_extr(:,3));
                    % Rescale! This is activity per pixel!               
                    activity(activity_i, seg_index,time_index) = av_activity / activity_num(activity_i, seg_index,time_index);
                else
                    activity(activity_i, seg_index,time_index)     = NaN;
                    activity_num(activity_i, seg_index,time_index) = 0;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%% SHIFTED ACTIVITY      %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if WINDOWS_TYPE == 1
                % Create mask
                img=zeros(imgSize);
                mask_image=roipoly(img,xb_shift,yb_shift);

                % Extract from the score map
                shift_activity_map_crop = double(mask_image) .* double(activity_map{activity_i});

                clear shift_activity_map_crop_extr;
                % Find the scores in the extracted segment of the score map
                [shift_activity_map_crop_extr(:,2), shift_activity_map_crop_extr(:,1), shift_activity_map_crop_extr(:,3)] = find(shift_activity_map_crop);

                % Get number of scores in segment
                shift_activity_num(activity_i, seg_index, time_index) = size(shift_activity_map_crop_extr,1);

                % Average shifted activity map
                if shift_activity_num(activity_i, seg_index, time_index) > 0
                    av_shift_activity = sum(shift_activity_map_crop_extr(:,3));
                    % Rescale! This is activity per pixel!
                    shift_activity(activity_i, seg_index,time_index) = av_shift_activity / shift_activity_num(activity_i, seg_index,time_index);

                else
                    shift_activity(activity_i, seg_index,time_index)     = NaN;
                    shift_activity_num(activity_i, seg_index,time_index) = 0;
                end
                end %WINDOWS_TYPE == 1
            end % DO_ACTIVITY
        end % for
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        

      
        
        %%%%%%%%%%%%%%%%%%  Nice plot    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PLOT_VECTORS
           relFrmNo = maskImgIndex(time)-firstImgIndex+1;
           if seg == START_SEG && seg_index == 1
              if ~exist('h_nice','var')
                 h_nice = figure;
              else
                 if ~ishandle(h_nice)
                    h_nice = figure;
                 end
              end

              figure(h_nice); hold off;
              %overlay_image = imread(overlay_image_path);
              %overlay_image = imread([IMG_DIR filesep filelist_images(time_index+prot_parameters.first_img-1,:)]);
              overlay_image = imread([IMG_DIR filesep filelist_images(relFrmNo,:)]);
              imshow(overlay_image,[0 1900]);
              hold on;
              titleStr = sprintf('Image Index: %d',maskImgIndex(time));
              title(titleStr);
              plot(x_p_edge, y_p_edge,'r');
              h_axes = gca;
              set(h_axes,'xgrid','off','ygrid','off','box','off','Visible','off');
              axis equal
              axis([0 imgSize(2) 0 imgSize(1)])
              axis ij
           end
           figure(h_nice);
           plot(xb,yb,'y','LineWidth',1);
           plot(xb_shift,  yb_shift,'y','LineWidth',1);
           if seg == START_SEG
              text(xb(1),yb(1),num2str(START_SEG),'Color','w');
           elseif seg == SEG_NR - END_SEG
              text(xb(end),yb(end), num2str(SEG_NR - END_SEG),'Color','w'); 
           end

           if DO_FLOW
              quiver('v6',av_velocity_field_crop(2), av_velocity_field_crop(1),...
                 av_velocity_field_crop(4), av_velocity_field_crop(3),20,'r');
              quiver('v6',av_shift_velocity_field_crop(2), av_shift_velocity_field_crop(1),...
                 av_shift_velocity_field_crop(4), av_shift_velocity_field_crop(3),20,'r');
           end

           if DO_SCORES
              plot(score_data_extr(:,1),score_data_extr(:,2),'rx');
           end

           if DO_PROT
              quiver(x_av_pos_prot(seg), y_av_pos_prot(seg), x_av_prot(seg), y_av_prot(seg), 10,'g');
           end
           %quiver(x_av_pos_normal(seg), y_av_pos_normal(seg), x_av_normal(seg), y_av_normal(seg),15,'m');

           %Modified by Lin Ji on Sept. 08, 2006.
           if seg == SEG_NR - END_SEG
              if ~isdir([merg_dir filesep 'overlayFig']);
                 success = mkdir(merg_dir,'overlayFig');
                 if ~success
                    error('Trouble making directory.');
                 end
              end
              hh = getframe(h_nice);
              %print(h_nice, [merg_dir filesep 'overlay' filesep 'overlay_' filelist_images(relFrmNo,:)],'-dtiff');
              imwrite(hh.cdata,[merg_dir filesep 'overlay' filesep 'overlay_' filelist_images(relFrmNo,:)],'tif');
              print(h_nice, [merg_dir 'overlay' filesep 'overlay_' filelist_images(relFrmNo,:) '.eps'],'-depsc2','-tiff');
              saveas(h_nice, [merg_dir filesep 'overlayFig' filesep 'overlay_' filelist_images(relFrmNo,1:end-3) 'fig'],'fig');
              %print(h_nice, [merg_dir 'nice.eps'],'-depsc2','-tiff');
              %hgsave(h_nice,[merg_dir 'nice.fig']);   
              %close(h_nice);
           end
        end % if nice plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    end % segment for
    
    % reset counter in matlab display
    fprintf(1,backSpc);
    
end %time  for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear score_matrix;
clear flow;
clear scores;
clear img;
clear mask_image;
clear mask_img;
clear score_map;
clear score_map_crop;
clear img_score_map_crop;
clear velocity_field;
clear velocity_field_crop;
clear x_p_edge;
clear y_p_edge;

%if PLOT_VECTORS
%    return;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Postprocessing  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  For postprocessing resolve the problem of missing data points in the 
%%%  activity maps. These points are designated with NaN. 
%%%  1. Determine position of NaN's
%%%  2. If user chose interpolate, than interpolate
%%%     else if user chose to replace NaN, then NaN's are replaced with 
%%%     if no interploation was choosen, than Nan are replaced with zeros
%%%     Images are not treated since image data is defined in every pixel
%%%     thus no missing data should occure.
%%%
%%%  3. Protrusion, flow, score, image data is aligned with linear
%%%     interpolation
%%%  4. Variables are saved to disk
%%%
%%%  5. Control measurments are calculated.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%               Data interpolation                   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INTERPOLATE
    num_segments   = size(protrusion_seg,1);
    num_time       = size(protrusion_seg,2);  
    
    % Get first the indices of the invalid scores. This is used 
    % later to clean the num_scores!
    if DO_SCORES
        index_scores = isnan(score);
        index_scores_plus = isnan(score_plus);
        index_scores_minus = isnan(score_minus);
        [x_scores,y_scores] = find(index_scores);
        [x_scores_plus,y_scores_plus] = find(index_scores_plus);
        [x_scores_plus,y_scores_minus] = find(index_scores_minus);
        
        
        if WINDOWS_TYPE == 1
            index_shift_scores = isnan(shift_score);
            index_shift_scores_plus = isnan(shift_score_plus);
            index_shift_scores_minus = isnan(shift_score_minus);
            [x_shift_scores,y_shift_scores] = find(index_shift_scores); 
            [x_shift_scores_plus,y_shift_scores_plus] = find(index_shift_scores_plus); 
            [x_shift_scores_minus,y_shift_scores_minus] = find(index_shift_scores_minus); 
        end
    end
    
    if DO_ACTIVITY_1 | DO_ACTIVITY_2
        if WINDOWS_TYPE == 1
            index_shift_activity = isnan(shift_activity);
            [x_shift_activity, y_shift_activity] = find(index_shift_activity);
        end
    end
    
    
    if DO_FLOW
        index_retro = isnan(retrograde_flow_normal);
        [x_retro, y_retro]=find(index_retro);
        if WINDOWS_TYPE == 1
            index_retro = isnan(shift_retrograde_flow_normal);
            [x_shift_retro, y_shift_retro]=find(index_retro);
        end
    end
    
    if INTERPOLATE == 1 
        %%%%% Interplolate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DO_SCORES
            score(:,2:num_time) =           prInterpolateActivity(score(:,2:end));
            score(:,1)     = 0;
            score_plus(:,2:num_time) =           prInterpolateActivity(score_plus(:,2:end));
            score_plus(:,1)     = 0;
            score_minus(:,2:num_time) =           prInterpolateActivity(score_minus(:,2:end));
            score_minus(:,1)     = 0;
            if WINDOWS_TYPE == 1
                shift_score(:,2:num_time) =     prInterpolateActivity(shift_score(:,2:end));
                shift_score(:,1)     = 0;
                shift_score_plus(:,2:num_time) =     prInterpolateActivity(shift_score_plus(:,2:end));
                shift_score_plus(:,1)     = 0;
                shift_score_minus(:,2:num_time) =     prInterpolateActivity(shift_score_minus(:,2:end));
                shift_score_minus(:,1)     = 0;
            end
        end
 
        if DO_FLOW
            retrograde_flow =               prInterpolateActivity(retrograde_flow);
            retrograde_flow_normal =        prInterpolateActivity(retrograde_flow_normal);
            if WINDOWS_TYPE == 1
                shift_retrograde_flow =         prInterpolateActivity(shift_retrograde_flow);
                shift_retrograde_flow_normal =  prInterpolateActivity(shift_retrograde_flow_normal);
            end
        end
    elseif INTERPOLATE == 2
        %%%%% Interplolate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DO_SCORES
            score(:,2:num_time) =           prDirInterpolateActivity(score(:,2:end));
            score(:,1)     = 0;
            score_plus(:,2:num_time) =           prInterpolateActivity(score_plus(:,2:end));
            score_plus(:,1)     = 0;
            score_minus(:,2:num_time) =           prInterpolateActivity(score_minus(:,2:end));
            score_minus(:,1)     = 0;
            
            if WINDOWS_TYPE == 1
                shift_score(:,2:num_time) =     prDirInterpolateActivity(shift_score(:,2:end));
                shift_score(:,1)     = 0;
                shift_score_plus(:,2:num_time) =     prInterpolateActivity(shift_score_plus(:,2:end));
                shift_score_plus(:,1)     = 0;
                shift_score_minus(:,2:num_time) =     prInterpolateActivity(shift_score_minus(:,2:end));
                shift_score_minus(:,1)     = 0;
            end
        end
%         if DO_ACTIVITY_1
%             shift_activity               =  prInterpolateActivity(shift_activity);   
%         end        
        if DO_FLOW
            retrograde_flow =               prDirInterpolateActivity(retrograde_flow);
            shift_retrograde_flow =         prDirInterpolateActivity(shift_retrograde_flow);
            if WINDOWS_TYPE == 1
                retrograde_flow_normal =        prDirInterpolateActivity(retrograde_flow_normal);
                shift_retrograde_flow_normal =  prDirInterpolateActivity(shift_retrograde_flow_normal); 
            end
        end
    end % method
end % interpolate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  End data interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Remove wrong values from the data sets and correct the %%%%
%%%%% number of samples!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INTERPOLATE
    % we set the number of values at the interpolated positions
    % to one
    if DO_SCORES
        score(:,1) = 0;
        score_plus(:,1) = 0;
        score_minus(:,1) = 0;
        for i=1:length(x_scores)
            score_num(x_scores(i),y_scores(i))  = 1;
        end
        for i=1:length(x_scores_plus)
            score_num(x_scores_plus(i),y_scores_plus(i))  = 1;
        end
        for i=1:length(x_scores_minus)
            score_num(x_scores_minus(i),y_scores_minus(i))  = 1;
        end
        if WINDOWS_TYPE == 1
            shift_score(:,1) = 0;
            for i=1:length(x_shift_scores)
                shift_score_num(x_shift_scores(i),y_shift_scores(i))  = 1;
            end
            shift_score_plus(:,1) = 0;
            for i=1:length(x_shift_scores_plus)
                shift_score_num_plus(x_shift_scores_plus(i),y_shift_scores_plus(i))  = 1;
            end
            shift_score_minus(:,1) = 0;
            for i=1:length(x_shift_scores_minus)
                shift_score_num_minus(x_shift_scores_minus(i),y_shift_scores_minus(i))  = 1;
            end
        end
    end
    if DO_ACTIVITY_1
        if WINDOWS_TYPE == 1
            for i=1:length(x_shift_activity)
                activity_num(x_shift_activity(i),y_shift_activity(i))  = 1;
            end
        end
    end
    if DO_FLOW
        for i=1:length(x_retro)
            retrograde_flow_num(x_retro(i),y_retro(i))  = 1;
        end

        if WINDOWS_TYPE == 1
            for i=1:length(x_shift_retro)
                shift_retrograde_flow_num(x_shift_retro(i),y_shift_retro(i))  = 1;
            end
        end
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if there was no interpolation just replace the NaN with zeros
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if REPLACE_NAN
        if DO_SCORES
            index_scores = isnan(score);
            [x_scores,y_scores] = find(index_scores);
            index_scores_plus = isnan(score_plus);
            [x_scores_plus,y_scores_plus] = find(index_scores_plus);
            index_scores_minus = isnan(score_minus);
            [x_scores_minus,y_scores_minus] = find(index_scores_minus);

            if WINDOWS_TYPE == 1
                index_shift_scores = isnan(shift_score);
                [x_shift_scores,y_shift_scores]=find(index_shift_scores);
                index_shift_scores_plus = isnan(shift_score_plus);
                [x_shift_scores_plus,y_shift_scores_plus]=find(index_shift_scores_plus);
                index_shift_scores_minus = isnan(shift_score_minus);
                [x_shift_scores_minus,y_shift_scores_minus]=find(index_shift_scores_minus);
            end
            for i=1:length(x_scores)
                score(x_scores(i),y_scores(i))      = 0;
                score_num(x_scores(i),y_scores(i))  = 0;
            end
            for i=1:length(x_scores_plus)
                score_plus(x_scores_plus(i),y_scores_plus(i))      = 0;
                score_num_plus(x_scores_plus(i),y_scores_plus(i))  = 0;
            end
            for i=1:length(x_scores_minus)
                score_minus(x_scores_minus(i),y_scores_minus(i))      = 0;
                score_num_minus(x_scores_minus(i),y_scores_minus(i))  = 0;
            end
            
            if WINDOWS_TYPE == 1
                for i=1:length(x_shift_scores)
                    shift_score(x_shift_scores(i),y_shift_scores(i))      = 0;
                    shift_score_num(x_shift_scores(i),y_shift_scores(i))  = 0;
                end
                for i=1:length(x_shift_scores_plus)
                    shift_score_plus(x_shift_scores_plus(i),y_shift_scores_plus(i))      = 0;
                    shift_score_num_plus(x_shift_scores_plus(i),y_shift_scores_plus(i))  = 0;
                end
                for i=1:length(x_shift_scores_minus)
                    shift_score_minus(x_shift_scores_minus(i),y_shift_scores_minus(i))      = 0;
                    shift_score_num_minus(x_shift_scores_minus(i),y_shift_scores_minus(i))  = 0;
                end
            end
        end

        if DO_FLOW
            index_retro = isnan(retrograde_flow_normal);
            [x_retro, y_retro]=find(index_retro);
            for i=1:length(x_retro)
                retrograde_flow(x_retro(i),y_retro(i))          = 0;
                retrograde_flow_normal(x_retro(i),y_retro(i))   = 0;
                retrograde_flow_num(x_retro(i),y_retro(i))      = 0;
            end

            if WINDOWS_TYPE == 1
                index_shift_retro = isnan(shift_retrograde_flow_normal);
                [x_shift_retro, y_shift_retro]=find(index_shift_retro);
                for i=1:length(x_shift_retro)
                    shift_retrograde_flow(x_shift_retro(i),  y_shift_retro(i))      = 0;
                    shift_retrograde_flow_normal(x_shift_retro(i),y_shift_retro(i)) = 0;
                    shift_retrograde_flow_num(x_shift_retro(i),y_shift_retro(i))    = 0;
                end
            end
        end % retro

        for activity_i = 1:2
            if DO_ACTIVITY(activity_i)
                index_activity = isnan(activity(activity_i,:,:));
                [x_activity, y_activity]=find(squeeze(index_activity));
                for i=1:length(x_activity)
                    activity(activity_i, x_activity(i), y_activity(i))   = 0;
                    activity_num(activity_i, x_activity(i), y_activity(i))   = 0;
                end

                if WINDOWS_TYPE == 1
                    index_shift_activity = isnan(shift_activity(activity_i,:,:));
                    [x_shift_activity, y_shift_activity]=find(squeeze(index_shift_activity));
                    for i=1:length(x_shift_activity)
                        shift_activity(activity_i,x_shift_activity(i),y_shift_activity(i))   = 0;
                        shift_activity_num(activity_i,x_shift_activity(i),y_shift_activity(i))   = 0;
                    end
                    clear index_shift_activity;
                    clear x_shift_activity, y_shift_activity;       
                end
                
                clear index_activity;
                clear x_activity, y_activity;
            end
        end
    end
end % if interpolate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  End data cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Correct for sampling at different time points  %%%%%%
%%%%%%%%%%%% This means: shift activity half time step %%%%%%%%%%%
%
%   time step    1     2     3     4       T-1     T
%   protrusion      1     2     3              T-1
%   activity     1     2     3     4               T
%   score        1     2     3     4       T-1
%   flow            1     2     3              T-1
if SHIFT_SCORES_HALF_TIME_STEP
    for activity_i = 1:2
        if DO_ACTIVITY(activity_i)        % activity(activity_nr, seg, time)
            for i =1:size(activity,3)-1
                activity(activity_i,:,i)       = (activity(activity_i,:,i)+activity(activity_i,:,i+1))./2;
                if WINDOWS_TYPE == 1
                    shift_activity(activity_i,:,i) = (shift_activity(activity_i,:,i) + shift_activity(activity_i,:,i+1))./2;
                end
            end
        end
    end

    if DO_SCORES
        % score(seg,time)
        for i =1:size(score,2)-1
            score(:,i) = (score(:,i)+score(:,i+1))./2;
            if WINDOWS_TYPE == 1
                shift_score(:,i) = (shift_score(:,i)+shift_score(:,i+1))./2;
            end
        end
        for i =1:size(score_plus,2)-1
            score_plus(:,i) = (score_plus(:,i)+score_plus(:,i+1))./2;
            if WINDOWS_TYPE == 1
                shift_score_plus(:,i) = (shift_score_plus(:,i)+shift_score_plus(:,i+1))./2;
            end
        end
        for i =1:size(score_minus,2)-1
            score_minus(:,i) = (score_minus(:,i)+score_minus(:,i+1))./2;
            if WINDOWS_TYPE == 1
                shift_score_minus(:,i) = (shift_score_minus(:,i)+shift_score_minus(:,i+1))./2;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate average for the selected time window length and time step.
% Then, save the variable to disk 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DO_PROT
    numColumns = size(protrusion_seg,2);
    if protTimeWinL <= 2 
       %When window length is set to 1, it means the timing of the protrusion is
       %exactly the corresponding frame number not the middle of the frame interval.
       protrusion_seg    = protrusion_seg(:,1:protTimeStep:end);
       protrusion_normal = protrusion_normal(:,1:protTimeStep:end);
       protrusion_angle  = protrusion_angle(:,1:protTimeStep:end);

       protStartTime = FIRST_TIME;
       protTimePoints = [1:protTimeStep:numColumns] + FIRST_TIME - 1;
       if protTimeWinL == 2 
          %Note: When window length is 2, it refers to one frame interval. The timing is
          %therefore the middle of the interval.
          protStartTime  = protStartTime + 0.5;
          protTimePoints = protTimePoints + 0.5;
       end
    else
       protStartTime  = (1+protTimeWinL)/2 + FIRST_TIME - 1;
       protTimePoints = [];

       protrusion_segTmp = protrusion_seg;
       protrusion_normalTmp = protrusion_normal;
       protrusion_angleTmp = protrusion_angle;
       iCol = 1;
       kk   = 1;
       while iCol+protTimeWinL-2 <= numColumns
          protrusion_seg(:,kk)    = mean(protrusion_segTmp(:,iCol:iCol+protTimeWinL-2),2);
          protrusion_normal(:,kk) = mean(protrusion_normalTmp(:,iCol:iCol+protTimeWinL-2),2);
          protrusion_angle(:,kk)    = mean(protrusion_angleTmp(:,iCol:iCol+protTimeWinL-2),2);

          protTimePoints = [protTimePoints iCol+(protTimeWinL-1)/2];

          iCol = iCol + protTimeStep;
          kk   = kk+1;
       end
       protrusion_seg(:,kk:end)    = [];
       protrusion_normal(:,kk:end) = [];
       protrusion_angle(:,kk:end)    = [];
    end

    save([merg_dir 'protrusion.mat'], 'protrusion_seg');
    save([merg_dir 'protrusion_n.mat'], 'protrusion_normal');
    
    displacement_seg = protrusion_seg.*TIME_INTERVAL;
    displacement_normal = protrusion_normal.*TIME_INTERVAL;
    
    save([merg_dir 'displacement.mat'], 'displacement_seg');
    save([merg_dir 'displacement_n.mat'], 'displacement_normal');

    merg_parameters.protStartTime  = protStartTime;
    merg_parameters.protTimePoints = protTimePoints;
end
if DO_SCORES
    numColumns = size(score,2);
    if scoreTimeWinL== 1
       %Always remove the first column so that the timing aligns with protrusion and flow
       % data. Since the poly/depoly score for each frame is calculated from 3
       % consecutive frames, the score of the very first frame is zero. By default,
       % we also calculate the average over two consecutive time intervals for the
       % flow and protrusion so that the activities synchronize in time. Otherwise,
       % there is a half time step mismatch.
       score = score(:,2:scoreTimeStep:end);
       score_plus = score_plus(:,2:scoreTimeStep:end);
       score_minus = score_minus(:,2:scoreTimeStep:end);
       if WINDOWS_TYPE == 1
          shift_score = shift_score(:,2:scoreTimeStep:end);
          shift_score_plus = shift_score_plus(:,2:scoreTimeStep:end);
          shift_score_minus = shift_score_minus(:,2:scoreTimeStep:end);
       end

       scoreStartTime  = FIRST_TIME+1;
       scoreTimePoints = [2:scoreTimeStep:numColumns] + FIRST_TIME - 1;
    else
       scoreStartTime  = (1+scoreTimeWinL)/2 + FIRST_TIME - 1;
       scoreTimePoints = []; 

       scoreTmp = score;
       scoreTmp_plus = score_plus;
       scoreTmp_minus = score_minus;
       if WINDOWS_TYPE == 1
          shift_scoreTmp = shift_score;
          shift_scoreTmp_plus = shift_score_plus;
          shift_scoreTmp_minus = shift_score_minus;
       end
       iCol = 1;
       kk   = 1;
       while iCol+scoreTimeWinL-1 <= numColumns
          score(:,kk)     = mean(scoreTmp(:,iCol:iCol+scoreTimeWinL-1),2);
          scoreTimePoints = [scoreTimePoints iCol+(scoreTimeWinL-1)/2];
          score_plus(:,kk)     = mean(scoreTmp_plus(:,iCol:iCol+scoreTimeWinL-1),2);
          score_minus(:,kk)     = mean(scoreTmp_minus(:,iCol:iCol+scoreTimeWinL-1),2);
          scoreTimePoints = [scoreTimePoints iCol+(scoreTimeWinL-1)/2];
          if WINDOWS_TYPE == 1
             shift_score(:,kk) = mean(shift_scoreTmp(:,iCol:iCol+scoreTimeWinL-1),2);
             shift_score_plus(:,kk) = mean(shift_scoreTmp_plus(:,iCol:iCol+scoreTimeWinL-1),2);
             shift_score_minus(:,kk) = mean(shift_scoreTmp_minus(:,iCol:iCol+scoreTimeWinL-1),2);
          end

          iCol = iCol + scoreTimeStep;
          kk   = kk+1;
       end
       score(:,kk:end)       = [];
       score_plus(:,kk:end)       = [];
       score_minus(:,kk:end)       = [];
       if WINDOWS_TYPE == 1
          shift_score(:,kk:end) = [];
          shift_score_plus(:,kk:end) = [];
          shift_score_minus(:,kk:end) = [];
       end
    end
    save([merg_dir 'turnover.mat'], 'score');
    save([merg_dir 'poly.mat'], 'score_plus');
    save([merg_dir 'depoly.mat'], 'score_minus');
    if WINDOWS_TYPE == 1
        save([merg_dir 'shift_turnover.mat'], 'shift_score');
        save([merg_dir 'shift_poly.mat'], 'shift_score_plus');
        save([merg_dir 'shift_depoly.mat'], 'shift_score_minus');
    end
    save([merg_dir 'turnover_num.mat'], 'score_num'); 
    save([merg_dir 'poly_num.mat'], 'score_num_plus');
    save([merg_dir 'depoly_num.mat'], 'score_num_minus');

    merg_parameters.scoreStartTime  = scoreStartTime;
    merg_parameters.scoreTimePoints = scoreTimePoints;
end
if DO_FLOW
    numColumns = size(retrograde_flow,2);
    if flowTimeWinL <= 2 
       %When window length is set to 1, it means the timing of the protrusion is
       %exactly the corresponding frame number not the middle of the frame interval.
       retrograde_flow        = retrograde_flow(:,1:flowTimeStep:end);
       retrograde_flow_normal = retrograde_flow_normal(:,1:flowTimeStep:end);
       retrograde_flow_angle  = retrograde_flow_angle(:,1:flowTimeStep:end);
       av_vector_field   = av_vector_field(:,1:flowTimeStep:end,:);
       if WINDOWS_TYPE == 1
          shift_retrograde_flow        = shift_retrograde_flow(:,1:flowTimeStep:end);
          shift_retrograde_flow_normal = shift_retrograde_flow_normal(:,1:flowTimeStep:end);
          av_shift_vector_field        = av_shift_vector_field(:,1:flowTimeStep:end,:);

          av_network_velocity_direction_cor = av_network_velocity_dirction_cor(:,1:flowTimeStep:end);
       end

       flowStartTime  = FIRST_TIME;
       flowTimePoints = [1:flowTimeStep:numColumns] + FIRST_TIME - 1;
       if flowTimeWinL == 2
          %Note: When window length is 2, it refers to one frame interval. The timing is
          %therefore the middle of the interval.
          flowStartTime  = flowStartTime + 0.5;
          flowTimePoints = flowTimePoints + 0.5;
       end
    else
       flowStartTime  = (1+flowTimeWinL)/2 + FIRST_TIME - 1;
       flowTimePoints = [];

       retrograde_flowTmp       = retrograde_flow;
       retrograde_flow_nTmp     = retrograde_flow_normal;
       retrograde_flow_angleTmp = retrograde_flow_angle;
       av_vector_fieldTmp       = av_vector_field;
       if WINDOWS_TYPE == 1
          shift_retrograde_flowTmp   = shift_retrograde_flow;
          shift_retrograde_flow_nTmp = shift_retrograde_flow_normal;
          av_shift_vector_fieldTmp   = av_shift_vector_field;

          av_network_velocity_direction_corTmp = av_network_velocity_direction_cor;
       end
       iCol = 1;
       kk   = 1;
       while iCol+flowTimeWinL-2 <= numColumns
          retrograde_flow(:,kk)   = mean(retrograde_flowTmp(:,iCol:iCol+flowTimeWinL-2),2);
          retrograde_flow_normal(:,kk) = mean(retrograde_flow_nTmp(:,iCol:iCol+flowTimeWinL-2),2);
          retrograde_flow_angle(:,kk)  = mean(retrograde_flow_angleTmp(:,iCol:iCol+flowTimeWinL-2),2);
          av_vector_field(:,kk,:) = mean(av_vector_fieldTmp(:,iCol:iCol+flowTimeWinL-2,:),2);
          if WINDOWS_TYPE == 1
             shift_retrograde_flow(:,kk)        = mean(shift_retrograde_flowTmp(:,iCol:iCol+flowTimeWinL-2),2);
             shift_retrograde_flow_normal(:,kk) = mean(shift_retrograde_flow_nTmp(:,iCol:iCol+flowTimeWinL-2),2);
             av_shift_vector_field(:,kk,:)      = mean(av_shift_vector_fieldTmp(:,iCol:iCol+flowTimeWinL-2,:),2);

             av_network_velocity_direction_cor(:,kk) = mean(av_network_velocity_direction_corTmp(:,iCol:iCol+flowTimeWinL-2),2);
          end

          flowTimePoints = [flowTimePoints iCol+(flowTimeWinL-1)/2];

          iCol = iCol + flowTimeStep;
          kk   = kk+1;
       end
       retrograde_flow(:,kk:end)        = [];
       retrograde_flow_normal(:,kk:end) = [];
       retrograde_flow_angle(:,kk:end)  = [];
       if WINDOWS_TYPE == 1
          shift_retrograde_flow(:,kk:end)        = [];
          shift_retrograde_flow_normal(:,kk:end) = [];
          av_shift_vector_field(:,kk:end,:)      = [];

          av_network_velocity_direction_cor(:,kk:end) = [];
       end
    end
    save([merg_dir 'av_vector_field'],'av_vector_field');
    save([merg_dir 'retrograde_flow'], 'retrograde_flow');
    save([merg_dir 'retrograde_flow_n'], 'retrograde_flow_normal');
    if WINDOWS_TYPE == 1
        save([merg_dir 'av_shift_vector_field'],'av_shift_vector_field');
        save([merg_dir 's_retrograde_flow'], 'shift_retrograde_flow');
      	save([merg_dir 's_retrograde_flow_n'], 'shift_retrograde_flow_normal'); 
        save([merg_dir 'av_network_velocity_direction_cor'], 'av_network_velocity_direction_cor');
    end

    merg_parameters.flowStartTime  = flowStartTime;
    merg_parameters.flowTimePoints = flowTimePoints;
end
if DO_ACTIVITY_1 | DO_ACTIVITY_2 
    numColumns = size(activity,3);
    if activityTimeWinL== 1
       %Always remove the first column so that the timing aligns with protrusion and flow
       % data. By default, we calculate the average over two consecutive time intervals for the
       % flow and protrusion so that the activities synchronize in time. Otherwise,
       % there is a half time step mismatch.
       activity = activity(:,:,2:activityTimeStep:end);
       if WINDOWS_TYPE == 1
          shift_activity = shift_activity(:,:,2:activityTimeStep:end);
       end

       activityStartTime  = FIRST_TIME+1;
       activityTimePoints = [2:activityTimeStep:numColumns] + FIRST_TIME - 1;
    else
       activityStartTime  = (1+activityTimeWinL)/2 + FIRST_TIME - 1;
       activityTimePoints = []; 

       activityTmp = activity;
       if WINDOWS_TYPE == 1
          shift_activityTmp = shift_activity;
       end
       iCol = 1;
       kk   = 1;
       while iCol+activityTimeWinL-1 <= numColumns
          activity(:,:,kk) = mean(activityTmp(:,:,iCol:iCol+activityTimeWinL-1),3);
          activityTimePoints = [activityTimePoints iCol+(activityTimeWinL-1)/2];
          if WINDOWS_TYPE == 1
             shift_activity(:,kk) = mean(shift_activityTmp(:,iCol:iCol+activityTimeWinL-1),2);
          end

          iCol = iCol + activityTimeStep;
          kk   = kk+1;
       end
       activity(:,:,kk:end)       = [];
       if WINDOWS_TYPE == 1
          shift_activity(:,:,kk:end) = [];
       end
    end
    if DO_ACTIVITY_1
        activity1_tmp = squeeze(activity(1,:,:));
        save([merg_dir 'activity1.mat'], 'activity1_tmp');
        if WINDOWS_TYPE == 1
            shift_activity1_tmp = squeeze(shift_activity(1,:,:));
            save([merg_dir 's_activity1.mat'], 'shift_activity1_tmp');
        end
    end
    
    if DO_ACTIVITY_2
        activity2_tmp = squeeze(activity(2,:,:));
        save([merg_dir 'activity2.mat'], 'activity2_tmp');
        if WINDOWS_TYPE == 1
            shift_activity2_tmp = squeeze(shift_activity(2,:,:));
            save([merg_dir 's_activity2.mat'], 'shift_activity2_tmp');
        end
    end

    merg_parameters.activityStartTime  = activityStartTime;
    merg_parameters.activityTimePoints = activityTimePoints;
end
%if DO_ACTIVITY_2
%    activity_tmp = squeeze(activity(2,:,:));
%    if WINDOWS_TYPE == 1
%        shift_activity_tmp = squeeze(shift_activity(2,:,:));
%    end
%    numColumns = size(activity_tmp,2);
%    if activityTimeWinL== 1
%       %Always remove the first column so that the timing aligns with protrusion and flow
%       % data. By default, we calculate the average over two consecutive time intervals for the
%       % flow and protrusion so that the activities synchronize in time. Otherwise,
%       % there is a half time step mismatch.
%       activity_tmp = activity_tmp(:,2:activityTimeStep:end);
%       if WINDOWS_TYPE == 1
%          shift_activity_tmp = shift_activity_tmp(:,2:activityTimeStep:end);
%       end
%
%       activity_2_StartTime  = FIRST_TIME+1;
%       activity_2_TimePoints = 2:activityTimeStep:numColumns + FIRST_TIME - 1;
%    else
%       activity_2_StartTime  = (1+activityTimeWinL)/2 + FIRST_TIME - 1;
%       activity_2_TimePoints = []; 
%
%       activity_tmpTmp = activity_tmp;
%       if WINDOWS_TYPE == 1
%          shift_activity_tmpTmp = shift_activity_tmp;
%       end
%       iCol = 1;
%       kk   = 1;
%       while iCol+activityTimeWinL-1 <= numColumns
%          activity_tmp(:,kk) = mean(activity_tmpTmp(:,iCol:iCol+activityTimeWinL-1),2);
%          activity_2_TimePoints = [activity_2_TimePoints iCol+(activityTimeWinL-1)/2];
%          if WINDOWS_TYPE == 1
%             shift_activity_tmp(:,kk) = mean(shift_activity_tmpTmp(:,iCol:iCol+activityTimeWinL-1),2);
%          end
%
%          iCol = iCol + activityTimeStep;
%          kk   = kk+1;
%       end
%       activity_tmp(:,kk:end)       = [];
%       if WINDOWS_TYPE == 1
%          shift_activity_tmp(:,kk:end) = [];
%       end
%    end
%    save([merg_dir 'activity2.mat'], 'activity_tmp');   
%    if WINDOWS_TYPE == 1
%        shift_activity_tmp = squeeze(shift_activity(2,:,:));    
%        save([merg_dir 's_activity2.mat'], 'shift_activity_tmp');   
%    end
%
%    merg_parameters.activity_2_StartTime  = activity_2_StartTime;
%    merg_parameters.activity_2_TimePoints = activity_2_TimePoints;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Write all the parameters to a file  %%%%%%%%%%%%%%%%%%%%%%%%%%%
save_pralpha_parameters([merg_dir 'merg_parameters.dat'], merg_parameters);
save([merg_dir  'merg_parameters'], 'merg_parameters');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Calculation of control values      %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   These does not affect the merged data, it is just used   %%
%%%%%%%%%%%%%   to double check the merging        %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% The variables have the structure [segments,time] %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Calculate TIME averaged values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Good for cell with strong variations in segments %%%%%%%%%%%%%%%%%%%
if DO_PROT
    av_protrusion_seg               = mean(protrusion_seg,2);
    av_protrusion_n_seg             = mean(protrusion_normal,2);
    av_protrusion_angle_seg         = mean(protrusion_angle,2);
    av_displacement_seg             = mean(displacement_seg,2);
    av_displacement_n_seg            = mean(displacement_normal,2);
    % calculate rms of protrusion
    for i=1:size(protrusion_normal,1)
        rms_protrusion_n_seg(i)    = norm(protrusion_normal(i,:))/sqrt(size(protrusion_normal,2));  
    end
end
if DO_FLOW  
    av_retrograde_flow_seg          = nanmean(retrograde_flow,2);
    av_retrograde_flow_n_seg        = nanmean(retrograde_flow_normal,2);
    av_num_seg_velocities_seg       = nanmean(retrograde_flow_num,2);
    if WINDOWS_TYPE == 1
       av_shift_retrograde_flow_n_seg  = nanmean(shift_retrograde_flow_normal,2);
       av_num_seg_shift_velocities_seg = nanmean(shift_retrograde_flow_num,2);
    end
    av_retrograde_flow_angle_seg    = nanmean(retrograde_flow_angle,2);
    
    av_vector_field_seg(:,1)        = nanmean(av_vector_field(:,1),2);
    av_vector_field_seg(:,2)        = nanmean(av_vector_field(:,2),2);
    if WINDOWS_TYPE == 1  
       av_shift_vector_field_seg(:,1)  = nanmean(av_shift_vector_field(:,1),2);
       av_shift_vector_field_seg(:,2)  = nanmean(av_shift_vector_field(:,2),2);
    end
end
if DO_SCORES
    av_turnover_seg                   = nanmean(score,2);
    av_turnover_sample_num_seg         = nanmean(score_num,2);
    av_poly_seg                   = nanmean(score_plus,2);
    av_poly_sample_num_seg         = nanmean(score_num_plus,2);
    av_depoly_seg                   = nanmean(score_minus,2);
    av_depoly_sample_num_seg         = nanmean(score_num_minus,2);
    total_turnover_sample_num_seg         = nansum(score_num,2);
    total_poly_sample_num_seg         = nansum(score_num_plus,2);
    total_depoly_sample_num_seg     = nansum(score_num_minus,2);
    
    if WINDOWS_TYPE == 1   
        av_shift_turnover_seg             = nanmean(shift_score,2);
        av_shift_turnover_sample_num_seg   = nanmean(shift_score_num,2);
        av_shift_poly_seg             = nanmean(shift_score_plus,2);
        av_shift_poly_sample_num_seg   = nanmean(shift_score_num_plus,2);
        av_shift_depoly_seg             = nanmean(shift_score_minus,2);
        av_shift_depoly_sample_num_seg   = nanmean(shift_score_num_minus,2);
        total_shift_turnover_sample_num_seg   = nansum(shift_score_num,2);
        total_shift_poly_sample_num_seg   = nansum(shift_score_num_plus,2);
        total_shift_depoly_sample_num_seg   = nansum(shift_score_num_minus,2);
        
        
    end
end
for activity_i = 1:2  
    if DO_ACTIVITY(activity_i)
    av_activity_seg(activity_i,:)                  = sum(activity(activity_i,:,:),3)     ./ sum(activity(activity_i,:,:) ~=0, 3);
    av_activity_sample_num_seg(activity_i,:)       = sum(activity_num(activity_i,:,:),3) ./ sum(activity_num(activity_i,:,:) ~=0, 3);    
    if WINDOWS_TYPE == 1 
    av_shift_activity_seg(activity_i,:)            = sum(shift_activity(activity_i,:,:),3)     ./ sum(shift_activity(activity_i,:,:) ~=0, 3);
    av_shift_activity_sample_num_seg(activity_i,:) = sum(shift_activity_num(activity_i,:,:),3) ./ sum(shift_activity_num(activity_i,:,:) ~=0, 3);     
    end
    end
end
% replace the Nan's
for activity_i = 1:2
    if DO_ACTIVITY(activity_i)
        index_activity = isnan(av_activity_seg(activity_i,:));
        [x_activity]=find(squeeze(index_activity));
        for i=1:length(x_activity)
            av_activity_seg(activity_i, x_activity(i))   = 0;
            av_activity_sample_num_seg(activity_i, x_activity(i))   = 0;
        end
        if WINDOWS_TYPE == 1
        index_shift_activity = isnan(av_shift_activity_seg(activity_i,:));
        [x_shift_activity]=find(squeeze(index_shift_activity));
        for i=1:length(x_shift_activity)
            av_shift_activity_seg(activity_i,x_shift_activity(i))   = 0;
            av_shift_activity_sample_num_seg(activity_i,x_shift_activity(i))   = 0;
        end
        clear index_shift_activity;
        clear x_shift_activity;        
        end
        clear index_activity;
        clear x_activity;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Calculate SEGMENT averaged values %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Good for cell with strong variations in time %%%%%%%%%%%%%%
if DO_PROT
    av_protrusion_time              = sum(protrusion_seg,1)      ./ sum(protrusion_seg ~=0, 1);
    av_protrusion_n_time            = sum(protrusion_normal,1)   ./ sum(protrusion_normal ~=0, 1);
    av_protrusion_angle_time        = sum(protrusion_angle,1)    ./ sum(protrusion_angle ~=0, 1);
    av_displacement_time            = sum(protrusion_seg*TIME_INTERVAL,1)    ./sum(protrusion_seg ~=0,1); 
    av_displacement_n_time          = sum(protrusion_normal*TIME_INTERVAL,1) ./sum(protrusion_normal ~=0,1);
end
if DO_FLOW
    av_retrograde_flow_time          = nanmean(retrograde_flow,1);
    av_retrograde_flow_n_time        = nanmean(retrograde_flow_normal,1);
    av_num_seg_velocities_time       = nanmean(retrograde_flow_num,1);
    if WINDOWS_TYPE == 1
       av_shift_retrograde_flow_n_time  = nanmean(shift_retrograde_flow_normal,1);       
       av_num_seg_shift_velocities_time = nanmean(shift_retrograde_flow_num,1);
    end
    av_retrograde_flow_angle_time    = nanmean(retrograde_flow_angle,1);
    
    av_vector_field_time(:,1)        = nanmean(av_vector_field(:,1),1);
    av_vector_field_time(:,1)        = nanmean(av_vector_field(:,2),1);
    if WINDOWS_TYPE == 1
       av_shift_vector_field_time(:,1)  = nanmean(av_shift_vector_field(:,1),1);
       av_shift_vector_field_time(:,1)  = nanmean(av_shift_vector_field(:,2),1);   
    end
end
if DO_SCORES
    av_turnover_time                   = nanmean(score,1);
    av_turnover_sample_num_time         = nanmean(score_num,1);
    av_poly_time                  = nanmean(score_plus,1);
    av_poly_sample_num_time         = nanmean(score_num_plus,1);
    av_depoly_time                  = nanmean(score_minus,1);
    av_depoly_sample_num_time         = nanmean(score_num_minus,1);
    total_turnover_sample_num_time         = nansum(score_num,1);
    total_poly_sample_num_time         = nansum(score_num_plus,1);
    total_depoly_sample_num_time         = nansum(score_num_minus,1);
    
    if WINDOWS_TYPE == 1
       av_shift_turnover_time             = nanmean(shift_score,1);
       av_shift_turnover_sample_num_time   = nanmean(shift_score_num,1);  
       av_shift_poly_time            = nanmean(shift_score_plus,1);
       av_shift_poly_sample_num_time  = nanmean(shift_score_num_plus,1);
       av_shift_depoly_time            = nanmean(shift_score_minus,1);
       av_shift_depoly_sample_num_time   = nanmean(shift_score_num_minus,1);
       total_shift_turnover_sample_num_time   = nansum(shift_score_num,1);  
       total_shift_poly_sample_num_time  = nansum(shift_score_num_plus,1);
       total_shift_depoly_sample_num_time   = nansum(shift_score_num_minus,1);
       
    end
end
for activity_i = 1:2  
    if DO_ACTIVITY(activity_i)
       av_activity_time(activity_i,:)                  = sum(activity(activity_i,:,:),2)      ./ sum(activity(activity_i,:,:) ~=0, 2);
       av_activity_sample_num_time(activity_i,:)       = sum(activity_num(activity_i,:,:),2)  ./ sum(activity_num(activity_i,:,:) ~=0, 2);
       if WINDOWS_TYPE == 1
          av_shift_activity_time(activity_i,:)            = sum(shift_activity(activity_i,:,:),2)     ./ sum(shift_activity(activity_i,:,:) ~=0, 2);
          av_shift_activity_sample_num_time(activity_i,:) = sum(shift_activity_num(activity_i,:,:),2) ./ sum(shift_activity_num(activity_i,:,:) ~=0, 2); 
       end
    end
end% replace the Nan's
for activity_i = 1:2
    if DO_ACTIVITY(activity_i)
        index_activity = isnan(av_activity_time(activity_i,:));
        [x_activity]=find(squeeze(index_activity));
        for i=1:length(x_activity)
            av_activity_time(activity_i, x_activity(i))   = 0;
            av_activity_sample_num_time(activity_i, x_activity(i))   = 0;
        end
        if WINDOWS_TYPE == 1
        index_shift_activity = isnan(av_shift_activity_time(activity_i,:));
        [x_shift_activity]=find(squeeze(index_shift_activity));
        for i=1:length(x_shift_activity)
            av_shift_activity_time(activity_i,x_shift_activity(i))   = 0;
            av_shift_activity_sample_num_time(activity_i,x_shift_activity(i))   = 0;
        end
        clear index_shift_activity;
        clear x_shift_activity;       
        end
        clear index_activity;
        clear x_activity;

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Calculate total averages (scalar values) %%%%%%%%%%%%%%%%%%%

% Note: av_total means average value of protrusion, displacement etc. per time
% and segement, i.e. average of the elements of protrusion, displacement
% etc. matrix over time and segment.

if DO_PROT
    av_total_prot_n             = mean(av_protrusion_n_seg);
    av_total_prot               = mean(av_protrusion_seg);
    av_total_disp_n             = mean(av_displacement_n_seg);
    av_total_disp               = mean(av_displacement_seg);
    av_net_normal_displacement             = mean(av_displacement_n_seg)*TOTAL_TIME;
    av_net_displacement               = mean(av_displacement_seg)*TOTAL_TIME;

end
if DO_FLOW
    av_total_retro_flow_n       = mean(av_retrograde_flow_n_seg);
    av_total_retro_flow         = mean(av_retrograde_flow_seg);
    if WINDOWS_TYPE == 1
    av_total_shift_retro_flow_n = mean(av_shift_retrograde_flow_n_seg);
    end
    av_total_vel_sample_num     = mean(av_num_seg_velocities_seg);   
end
if DO_SCORES
    av_total_turnover             = mean(av_turnover_seg );
    av_total_turnover_sample_num   = mean(av_turnover_sample_num_seg);
    av_total_poly             = mean(av_poly_seg );
    av_total_poly_sample_num   = mean(av_poly_sample_num_seg);
    av_total_depoly            = mean(av_depoly_seg );
    av_total_depoly_sample_num   = mean(av_depoly_sample_num_seg);
    
    if WINDOWS_TYPE == 1
    av_total_shift_turnover             = mean(av_shift_turnover_seg );
    av_total_shift_turnover_sample_num   = mean(av_shift_turnover_sample_num_seg);  
    av_total_shift_poly            = mean(av_shift_poly_seg );
    av_total_shift_poly_sample_num   = mean(av_shift_poly_sample_num_seg);
    av_total_shift_depoly           = mean(av_shift_depoly_seg );
    av_total_shift_depoly_sample_num   = mean(av_shift_depoly_sample_num_seg);
    end
end
for activity_i = 1:2
    if DO_ACTIVITY(activity_i)
    av_total_activity(activity_i)             = mean(av_activity_seg(activity_i,:));
    av_total_activity_sample_num(activity_i)  = mean(av_activity_sample_num_seg(activity_i,:)); 
    if WINDOWS_TYPE == 1
    av_total_shift_activity(activity_i)             = mean(av_shift_activity_seg(activity_i,:));
    av_total_shift_activity_sample_num(activity_i)  = mean(av_shift_activity_sample_num_seg(activity_i,:));
    end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the x-axis used for ploting the data
x_seg_axis  = START_SEG:SEG_NR - END_SEG;
%x_time_axis = (FIRST_TIME: 1: TOTAL_TIME_STEPS).* TIME_INTERVAL;


if DO_PROT 
   prot_x_time_axis = protTimePoints.* TIME_INTERVAL;
   x_time_axis = protTimePoints.* TIME_INTERVAL;
end
if DO_FLOW
   flow_x_time_axis = flowTimePoints.* TIME_INTERVAL;
   x_time_axis = flowTimePoints.* TIME_INTERVAL;
end
if DO_SCORES
   score_x_time_axis = scoreTimePoints.* TIME_INTERVAL;
   x_time_axis = scoreTimePoints.* TIME_INTERVAL;
end
if DO_ACTIVITY(1)
    act_x_time_axis{1} = activityTimePoints.* TIME_INTERVAL;
    x_time_axis = activityTimePoints.* TIME_INTERVAL;
end
if DO_ACTIVITY(2)
    act_x_time_axis{2} = activityTimePoints.* TIME_INTERVAL;
    x_time_axis = activityTimePoints.* TIME_INTERVAL;
end
if (TOTAL_TIME_STEPS - FIRST_TIME) > 5
    if DO_PROT
        %%%%%%%%%%%%%%   Calculate a robust estimage of protrusion %%%%%%%
        [robust_av_protrusion, stats] = robustfit(prot_x_time_axis, av_protrusion_n_time);
        robust_sigma_protrusion = stats.robust_s;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    if DO_FLOW
        %%%%%%%%%%%%%%   Calculate a robust estimage of retrograde flow %%
        [robust_av_retrograde_flow, stats] = robustfit(flow_x_time_axis, av_retrograde_flow_n_time);
        robust_sigma_retrograde_flow = stats.robust_s;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%   Calculate a robust estimage of shifted retrograde flow %%
        if WINDOWS_TYPE == 1
            [robust_av_shift_retrograde_flow, stats] = robustfit(flow_x_time_axis, av_shift_retrograde_flow_n_time);
            robust_sigma_shift_retrograde_flow = stats.robust_s;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    if DO_SCORES
        %%%%%%%%%%%%%%   Calculate a robust estimage of scores %%%%%%%%%%%
        [robust_av_turnover, stats] = robustfit(score_x_time_axis, av_turnover_time);
        robust_sigma_turnover = stats.robust_s; 
        [robust_av_poly, stats] = robustfit(score_x_time_axis, av_poly_time);
        robust_sigma_poly = stats.robust_s;
        [robust_av_depoly, stats] = robustfit(score_x_time_axis, av_depoly_time);
        robust_sigma_depoly = stats.robust_s;
        
        if WINDOWS_TYPE == 1
            [robust_av_shift_turnover, stats] = robustfit(score_x_time_axis, av_shift_turnover_time);
            robust_sigma_shift_turnover = stats.robust_s;
            [robust_av_shift_poly, stats] = robustfit(score_x_time_axis, av_shift_poly_time);
            robust_sigma_shift_poly = stats.robust_s;
            [robust_av_shift_depoly, stats] = robustfit(score_x_time_axis, av_shift_depoly_time);
            robust_sigma_shift_depoly = stats.robust_s;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    for activity_i = 1:2
       if DO_ACTIVITY(activity_i)
          %%%%%%%%%%%%%%   Calculate a robust estimage of activity%%%%%%%%%%%
          [robust_av_activity{activity_i}, stats(activity_i)] = robustfit(act_x_time_axis{activity_i}, av_activity_time(activity_i,:));
          robust_sigma_activity(activity_i) = stats(activity_i).robust_s; 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if WINDOWS_TYPE == 1
             [robust_av_shift_activity{activity_i}, stats(activity_i)] = robustfit(act_x_time_axis{activity_i}, av_shift_activity_time(activity_i,:));
             robust_sigma_shift_activity(activity_i) = stats(activity_i).robust_s;   
          end
       end
    end
else
    %no robust regression can be performed if there are too few
    %parameters!
    if DO_PROT
        robust_av_protrusion(1) =               av_total_prot_n;
        robust_av_protrusion(2) =               0;
        robust_sigma_protrusion =               0; 
    end
    if DO_FLOW
        robust_av_retrograde_flow(1) =          av_total_retro_flow_n;
        robust_av_retrograde_flow(2) =          0;
        robust_sigma_retrograde_flow =          0;

        if WINDOWS_TYPE == 1
        robust_av_shift_retrograde_flow(1) =    av_total_shift_retro_flow_n;
        robust_av_shift_retrograde_flow(2) =    0;
        robust_sigma_shift_retrograde_flow =    0;
        end
    end
    if DO_SCORES
        robust_av_turnover(1) =                av_total_turnover;
        robust_av_turnover(2) =                av_total_turnover;
        robust_sigma_turnover =                0;
        
        robust_av_poly(1) =                    av_total_poly;
        robust_av_poly(2) =                    av_total_poly;
        robust_sigma_poly =                    0;
        
        robust_av_depoly(1) =                  av_total_depoly;
        robust_av_depoly(2) =                  av_total_deply;
        robust_sigma_deploly =                 0;
        
        if WINDOWS_TYPE == 1
        robust_av_shift_turnover(1) =          av_total_shift_turnover;
        robust_av_shift_turnover(2) =          av_total_shift_turnover;
        robust_sigma_shift_turnover =          0;
        
        robust_av_shift_poly(1) =              av_total_shift_poly;
        robust_av_shift_poly(2) =              av_total_shift_poly;
        robust_sigma_shift_poly =              0;
        
        robust_av_shift_depoly(1) =            av_total_shift_depoly;
        robust_av_shift_depoly(2) =            av_total_shift_depoly;
        robust_sigma_shift_depoly =            0;
        
        
        end
    end    
    for activity_i = 1:2
    if DO_ACTIVITY(activity_i)
        robust_av_activity{activity_i} =        [av_total_activity;av_total_activity];
        robust_sigma_activity =                 0;
        if WINDOWS_TYPE == 1
        robust_av_shift_activity{activity_i} =  [av_total_shift_activity;av_total_shift_activity];
        robust_sigma_shift_activity =           0;        
        end
    end
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Save the average variables to disk %%%%%%%%%%
if DO_PROT
    save([av_data_dir 'av_prot_variables.mat'], 'av_protrusion_n_seg','av_protrusion_n_time',...
        'robust_av_protrusion', 'robust_sigma_protrusion');
    
    save([av_data_dir 'av_disp_variables.mat'], 'av_displacement_n_seg','av_displacement_n_time',...
        'av_total_disp_n', 'av_total_disp', 'av_net_normal_displacement', 'av_net_displacement');
    
end
if DO_FLOW
    save([av_data_dir 'av_retro_variables.mat'], 'av_retrograde_flow_n_seg', 'av_num_seg_velocities_seg',...
        'av_retrograde_flow_n_time','av_num_seg_velocities_time',...
        'robust_av_retrograde_flow', 'robust_sigma_retrograde_flow');

    if WINDOWS_TYPE == 1
    save([av_data_dir 'av_shift_retro_variables.mat'], 'av_shift_retrograde_flow_n_seg', 'av_num_seg_shift_velocities_time',...
        'av_shift_retrograde_flow_n_time','av_num_seg_shift_velocities_seg',...
        'robust_av_shift_retrograde_flow', 'robust_sigma_shift_retrograde_flow');
    end
end
if DO_SCORES
    save([av_data_dir 'av_turnover.mat'], 'av_turnover_seg', 'av_turnover_sample_num_seg',...
        'total_turnover_sample_num_seg','av_turnover_time','av_turnover_sample_num_time',...
        'total_turnover_sample_num_time','robust_av_turnover', 'robust_sigma_turnover');
    
    save([av_data_dir 'av_poly.mat'], 'av_poly_seg', 'av_poly_sample_num_seg',...
        'total_poly_sample_num_seg','av_poly_time','av_poly_sample_num_time',...
        'total_poly_sample_num_time','robust_av_poly', 'robust_sigma_poly');
    
    save([av_data_dir 'av_depoly.mat'], 'av_depoly_seg', 'av_depoly_sample_num_seg',...
        'total_depoly_sample_num_seg','av_depoly_time','av_depoly_sample_num_time',...
        'total_depoly_sample_num_time','robust_av_depoly', 'robust_sigma_depoly');
    
    
    if WINDOWS_TYPE == 1
    save([av_data_dir 'av_shift_turnover.mat'], 'av_shift_turnover_seg', 'av_shift_turnover_sample_num_seg',...
        'total_shift_turnover_sample_num_seg','av_shift_turnover_time','av_shift_turnover_sample_num_time',...
        'total_shift_turnover_sample_num_time','robust_av_shift_turnover', 'robust_sigma_shift_turnover');
    
    save([av_data_dir 'av_shift_poly.mat'], 'av_shift_poly_seg', 'av_shift_poly_sample_num_seg',...
        'total_shift_poly_sample_num_seg','av_shift_poly_time','av_shift_poly_sample_num_time',...
        'total_shift_poly_sample_num_time','robust_av_shift_poly', 'robust_sigma_shift_poly');
    
    save([av_data_dir 'av_shift_depoly.mat'], 'av_shift_depoly_seg', 'av_shift_depoly_sample_num_seg',...
        'total_shift_depoly_sample_num_seg','av_shift_depoly_time','av_shift_depoly_sample_num_time',...
        'total_shift_depoly_sample_num_time','robust_av_shift_depoly', 'robust_sigma_shift_depoly');
    
    
    end
end
if DO_ACTIVITY_1 || DO_ACTIVITY_2
    save([av_data_dir 'av_activity.mat'], 'av_activity_seg', 'av_activity_sample_num_seg',...
        'av_activity_time','av_activity_sample_num_time',...
        'robust_av_activity', 'robust_sigma_activity');
    if WINDOWS_TYPE == 1
    save([av_data_dir 'av_shift_activity.mat'], 'av_shift_activity_seg', 'av_shift_activity_sample_num_seg',...
        'av_shift_activity_time','av_shift_activity_sample_num_time',...
        'robust_av_shift_activity', 'robust_sigma_shift_activity');  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot TIME average values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if UNITS == 0
    xLabelText = 'frame #';
    yLabelText = 'Velocity (pixel/frame)';
    yLabelTextDisplacement = 'Displacement (pixel)';
    yLabelTextSegmentLength = 'Segment Length (pixel)';
elseif UNITS == 1
    xLabelText = 'Time (s)';    
    yLabelText = 'Velocity (nm/s)';
    yLabelTextDisplacement = 'Displacement (nm)';
    yLabelTextSegmentLength = 'Segment Length (nm)';
else
    xLabelText = 'Time (min)';    
    yLabelText = 'Velocity (um/min)';
    yLabelTextDisplacement = 'Displacement (Micron)';
    yLabelTextSegmentLength = 'Segment Length (Micron)'
end
%
% 0: sum scores in window
% 1: sum scores in window and convert it to /sec
% 2: average scores in window
% 3: average scores in window and convert it to /sec
if SCORES_CONVERT == 0 
    yScoreLabelText = 'Turnover: (sum of scores in window) ';
    yScoreLabelTextplus = 'Poly: (sum of poly in window) ';
    yScoreLabelTextminus = 'Depoly: (sum of depoly in window) ';
    
elseif SCORES_CONVERT == 1 
    yScoreLabelText = 'Turnover: (sum of scores in window / sec )';
    yScoreLabelTextplus = 'Poly: (sum of poly in window / sec )';
    yScoreLabelTextminus = 'Depoly: (sum of depoly in window / sec )';
    
elseif SCORES_CONVERT == 2 
    yScoreLabelText = 'Turnover: (ave. scores in window )';
    yScoreLabelTextplus = 'Polymerization: (ave. poly in window )';
    yScoreLabelTextminus = 'Polymerization: (ave. depoly in window )';
    
else
    yScoreLabelText = 'Turnover: (ave. scores in window / sec )'; 
    yScoreLabelTextplus = 'Polymerization: (ave. poly in window / sec )';
    yScoreLabelTextminus = 'Polymerization: (ave. depoly in window / sec )';
end



% Re-scale segment length
segment_length = segment_length .* PIXEL;
% Determine average segment length
segment_length_av = mean(segment_length);
figure
x_axis = (FIRST_TIME: 1: TOTAL_TIME_STEPS).* TIME_INTERVAL;
plot(x_axis, segment_length)
title('Average Segment Length');
xlabel(xLabelText);
ylabel(yLabelTextSegmentLength);
save([merg_dir 'segment_length_av'],'segment_length_av');




if WINDOWS_TYPE ~= 2  
    
if DO_FLOW && DO_SCORES && DO_PROT
    h_time_av1 = figure;
%     subplot(2,1,1);
%     plot(x_seg_axis, av_protrusion_n_seg,'g');
%     hold on
%     plot(x_seg_axis, av_retrograde_flow_n_seg,'r');  
%     plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_prot_n,          '-g');
%     text(START_SEG+1, 1.3*av_total_prot_n,        num2str(av_total_prot_n));
%     plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_retro_flow_n,    '-r');
%     text(START_SEG+1, 1.3*av_total_retro_flow_n,  num2str(av_total_retro_flow_n));
    
%     ylabel(yLabelText);
%     legend('Protrusion','Retrograde flow','Location','Best');
% 
%     subplot(2,1,2);
    plot(x_seg_axis, av_turnover_seg , 'b');
    hold on
    plot(x_seg_axis, av_shift_turnover_seg , '--b');
    plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_turnover,  'b');
    text(START_SEG+1, 1.4*av_total_turnover,        num2str(av_total_turnover));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yScoreLabelText);
    legend('Turnover lamellipodium', 'Turnover lamella', 'Location','Best');
    
    hgsave(h_time_av1,[figures_dir 'time_av_turnover.fig']);
    print(h_time_av1, [figures_dir 'time_av_turnover.eps'],'-depsc2','-tiff');
    print(h_time_av1, [figures_dir 'time_av_turnover.tif'],'-dtiff');
    
    h_time_av2 = figure;
    subplot(2,1,1);
    plot(x_seg_axis, av_poly_seg , 'b');
    hold on
    plot(x_seg_axis, av_shift_poly_seg , '--b');
    plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_poly,  'b');
    text(START_SEG+1, 1.4*av_total_poly,        num2str(av_total_poly));
    
    ylabel(yScoreLabelTextplus);
    legend('Poly in lamellipodium', 'Poly in lamella','Location','Best');
    title('Time statistics');
    
    subplot(2,1,2);
    plot(x_seg_axis, av_depoly_seg , 'b');
    hold on
    plot(x_seg_axis, av_shift_depoly_seg , '--b');
    plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_depoly,  'b');
    text(START_SEG+1, 1.4*av_total_depoly,        num2str(av_total_depoly));
    xlabel('Segment #');
    ylabel(yScoreLabelTextminus);
    legend('De-poly in lamellipodium', 'De-poly in lamella','Location','Best');
    
    hgsave(h_time_av2,[figures_dir 'time_av_poly-depoly.fig']);
    print(h_time_av2, [figures_dir 'time_av_poly-depoly.eps'],'-depsc2','-tiff');
    print(h_time_av2, [figures_dir 'time_av_poly-depoly.tif'],'-dtiff');
end

% Plot only the protrusion and retrograde flow 
if DO_PROT && ~DO_FLOW
    h_red_time_av = figure;
    plot(x_seg_axis, av_protrusion_n_seg,'g');
    hold on
    plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_prot_n,'-g');
    text(min(x_seg_axis), 1.4*av_total_prot_n,             num2str(av_total_prot_n));
    xlim([x_seg_axis(1) x_seg_axis(end)]);
    legend('Protrusion','Location','Best');
    title('Time statistics');
    xlabel('Segment #'); 
    ylabel(yLabelText);
    
    
    h_rms_prot = figure;
    plot(x_seg_axis, rms_protrusion_n_seg,'g');
    legend('RMS Protrusion', 'Location','Best');
    title('Time statistics');
    xlabel('Segment #'); 
    ylabel(yLabelText);
    
    
elseif DO_PROT && DO_FLOW
    h_red_time_av = figure;
    plot(x_seg_axis, av_protrusion_n_seg,'g');
    hold on
    plot(x_seg_axis, av_retrograde_flow_n_seg,'r');
    plot(x_seg_axis, av_shift_retrograde_flow_n_seg,  ':r');  
    legend('Protrusion','Retrograde flow','Shifted retrograde flow','Location','Best');
    text(min(x_seg_axis), 1.4*av_total_prot_n,             num2str(av_total_prot_n));
    plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_prot_n,'-g');    
    plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_retro_flow_n,'-r');
    plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_shift_retro_flow_n,':r');
    text(min(x_seg_axis), 1.4*av_total_retro_flow_n,       num2str(av_total_retro_flow_n));
    text(min(x_seg_axis), 1.4*av_total_shift_retro_flow_n, num2str(av_total_shift_retro_flow_n));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yLabelText);
end


if DO_PROT && DO_SCORES
    % Set the second y axis for scores
%     ax1 = gca;
%     set(ax1,'box','off');
%     xlimits1 = get(ax1,'XLim');
%     ylimits1 = get(ax1,'YLim');
%     ylimits2(1) =     min(av_turnover_seg);
%     ylimits2(2) =     max(av_turnover_seg);
%     if abs(ylimits1(1) / ylimits1(2)) < abs(ylimits2(1) / ylimits2(2))
%         ylimits2(2) = (ylimits1(2) / ylimits1(1)) * ylimits2(1);
%     else
%         ylimits2(1) = (ylimits1(1) / ylimits1(2)) * ylimits2(2);
%     end
%     ax2 = axes('position',get(ax1,'position'));
%     set(ax2,'XAxisLocation','bottom','YAxisLocation','right','color','no');
%     axis(ax2,[xlimits1(1) xlimits1(2) ylimits2(1)  ylimits2(2)]);

    % plot average scores
%     line(x_seg_axis ,av_turnover_seg ,'Color','b');
%     line(x_seg_axis ,av_shift_turnover_seg ,'Color','b','LineStyle','--');
%     line(x_seg_axis,  ones(size(x_seg_axis)).*av_total_turnover, 'Color', 'b');
%     text(min(x_seg_axis), 1.4*av_total_turnover,        num2str(av_total_turnover));
%     ylabel(yScoreLabelText);
    
    figure
    line(x_seg_axis ,av_turnover_seg ,'Color','b');
    line(x_seg_axis ,av_shift_turnover_seg ,'Color','b','LineStyle','--');
    line(x_seg_axis,  ones(size(x_seg_axis)).*av_total_turnover, 'Color', 'b');
    text(min(x_seg_axis), 1.4*av_total_turnover,        num2str(av_total_turnover));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yScoreLabelText);
    legend('Turnover in Lamellipodium','Turnover in Lamella', 'Location','Best');
    
    figure
    line(x_seg_axis ,av_poly_seg ,'Color','b');
    line(x_seg_axis ,av_shift_poly_seg ,'Color','b','LineStyle','--');
    line(x_seg_axis,  ones(size(x_seg_axis)).*av_total_poly, 'Color', 'b');
    text(min(x_seg_axis), 1.4*av_total_poly,        num2str(av_total_poly));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yScoreLabelTextplus);
    legend('Poly in Lamellipodium','Poly in Lamella', 'Location','Best');
    
    figure
    line(x_seg_axis ,av_depoly_seg ,'Color','b');
    line(x_seg_axis ,av_shift_depoly_seg ,'Color','b','LineStyle','--');
    line(x_seg_axis,  ones(size(x_seg_axis)).*av_total_depoly, 'Color', 'b');
    text(min(x_seg_axis), 1.4*av_total_depoly,        num2str(av_total_depoly));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yScoreLabelTextminus);
    legend('De-Poly in Lamellipodium','De-Poly in Lamella', 'Location','Best');
end

if DO_PROT
    hgsave(h_red_time_av,[figures_dir 'time_av_prot_retro.fig']); 
    print(h_red_time_av, [figures_dir 'time_av_prot_retro.eps'],'-depsc2','-tiff'); 
    print(h_red_time_av, [figures_dir 'time_av_prot_retro.tif'],'-dtiff'); 
end

if DO_SCORES && ~DO_PROT
    figure
    line(x_seg_axis ,av_turnover_seg ,'Color','b');
    line(x_seg_axis ,av_shift_turnover_seg ,'Color','b','LineStyle','--');
    line(x_seg_axis,  ones(size(x_seg_axis)).*av_total_turnover, 'Color', 'b');
    text(min(x_seg_axis), 1.4*av_total_turnover,        num2str(av_total_turnover));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yScoreLabelText);
    legend('Turnover in Lamellipodium','Turnover in Lamella', 'Location','Best');
    
    figure
    line(x_seg_axis ,av_poly_seg ,'Color','b');
    line(x_seg_axis ,av_shift_poly_seg ,'Color','b','LineStyle','--');
    line(x_seg_axis,  ones(size(x_seg_axis)).*av_total_poly, 'Color', 'b');
    text(min(x_seg_axis), 1.4*av_total_poly,        num2str(av_total_poly));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yScoreLabelTextplus);
    legend('Poly in Lamellipodium','Poly in Lamella', 'Location','Best');
    
    figure
    line(x_seg_axis ,av_depoly_seg ,'Color','b');
    line(x_seg_axis ,av_shift_depoly_seg ,'Color','b','LineStyle','--');
    line(x_seg_axis,  ones(size(x_seg_axis)).*av_total_depoly, 'Color', 'b');
    text(min(x_seg_axis), 1.4*av_total_depoly,        num2str(av_total_depoly));
    title('Time statistics');
    xlabel('Segment #');
    ylabel(yScoreLabelTextminus);
    legend('De-Poly in Lamellipodium','De-Poly in Lamella', 'Location','Best');
    
end

if DO_ACTIVITY_1
    h_activity_1_time_av = figure;
    plot(x_seg_axis, av_activity_seg(1,:),'m');
    hold on
    plot(x_seg_axis, av_shift_activity_seg(1,:),'m--');
    legend('Activity 1','Activity in the Lamellum 1','Location','Best');
    
    % plot average activities
    plot(x_seg_axis,  ones(size(x_seg_axis)).*av_total_activity(1), 'm');
    text(START_SEG+1, av_total_activity(1),        num2str(av_total_activity(1)));
    plot(x_seg_axis,  ones(size(x_seg_axis)).*av_total_shift_activity(1),  'm--');
    text(START_SEG+1, av_total_shift_activity(1),        num2str(av_total_shift_activity(1)));
    xlim([x_seg_axis(1) x_seg_axis(end)]);
    
    title('Time statistics 1');
    ylabel('Intensity per pixel');   
    xlabel('Segment #');   
    
    hgsave(h_activity_1_time_av,[figures_dir 'time_av_activity_1.fig']); 
    print(h_activity_1_time_av, [figures_dir 'time_av_activity_1.eps'],'-depsc2','-tiff'); 
    print(h_activity_1_time_av, [figures_dir 'time_av_activity_1.tif'],'-dtiff');
end
if DO_ACTIVITY_2
    h_activity_2_time_av = figure;
    plot(x_seg_axis, av_activity_seg(2,:),'m');
    hold on
    plot(x_seg_axis, av_shift_activity_seg(2,:),'m--');
    legend('Activity 2','Activity in the Lamellum 2','Location','Best');
    
    % plot average activities
    plot(x_seg_axis,  ones(size(x_seg_axis)).*av_total_activity(2), 'm');
    text(START_SEG+1, av_total_activity(2),        num2str(av_total_activity(2)));
    plot(x_seg_axis,  ones(size(x_seg_axis)).*av_total_shift_activity(2),  'm--');
    text(START_SEG+1, av_total_shift_activity(2),        num2str(av_total_shift_activity(2)));
    xlim([x_seg_axis(1) x_seg_axis(end)]);
    
    title('Time statistics 2');
    ylabel('Intensity per pixel');   
    xlabel('Segment #');   
    
    hgsave(h_activity_2_time_av,[figures_dir 'time_av_activity_2.fig']); 
    print(h_activity_2_time_av, [figures_dir 'time_av_activity_2.eps'],'-depsc2','-tiff'); 
    print(h_activity_2_time_av, [figures_dir 'time_av_activity_2.tif'],'-dtiff');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot the angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DO_PROT
    h_angle = figure;
end
if DO_PROT && DO_FLOW
    plot(x_seg_axis, av_retrograde_flow_angle_seg, '-sr')
    hold on
    plot(x_seg_axis, av_protrusion_angle_seg, '-db')
    xlabel('Segment #');
    ylabel('Angle [degree]');
    title('Time Statistics');
    legend('Retrograde flow - normal direction', 'Protrusion - normal direction', 'Location', 'Best');
elseif DO_PROT
    plot(x_seg_axis, av_protrusion_angle_seg, '-db')
    xlabel('Segment #');
    ylabel('Angle [degree]');
    title('Time statistics');
    legend( 'Protrusion - normal direction', 'Location', 'Best');
end
if DO_PROT
    hgsave(h_angle,[figures_dir 'angle.fig']); 
    print(h_angle, [figures_dir 'angle.eps'],'-depsc2','-tiff'); 
    print(h_angle, [figures_dir 'angle.tif'],'-dtiff');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % if WINDOWS_TYPE ~ 2  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot SEGMENT average values %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DO_SCORES && DO_FLOW
    h_seg_av1 = figure;
%     subplot(2,1,1);
%     plot(prot_x_time_axis, av_protrusion_n_time,'g');
%     hold on
%     plot(flow_x_time_axis, av_retrograde_flow_n_time,       'r');
%     plot(flow_x_time_axis, av_shift_retrograde_flow_n_time,  ':r');
%     plot(prot_x_time_axis, ones(size(flow_x_time_axis)).*av_total_prot_n,'-g');
%     legend('Protrusion','Retrograde flow','Shifted retrograde flow', 'Location', 'Best');
%     text(min(prot_x_time_axis), 1.4*av_total_prot_n,             num2str(av_total_prot_n));
% 
%     plot(flow_x_time_axis, av_retrograde_flow_n_time,       'r');
%     plot(flow_x_time_axis, av_shift_retrograde_flow_n_time,  ':r');
%     plot(flow_x_time_axis, ones(size(flow_x_time_axis)).*av_total_retro_flow_n,            '-r');
%     plot(flow_x_time_axis, ones(size(flow_x_time_axis)).*av_total_shift_retro_flow_n,      ':r');
%     text(min(flow_x_time_axis), 1.4*av_total_retro_flow_n,       num2str(av_total_retro_flow_n));
%     text(min(flow_x_time_axis), 1.4*av_total_shift_retro_flow_n, num2str(av_total_shift_retro_flow_n));

    
%     ylabel(yLabelText);
    
%     subplot(2,1,2);
    plot(score_x_time_axis,  av_turnover_time , 'b');
    hold on
    plot(score_x_time_axis,  av_shift_turnover_time , '--b');
    plot(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_turnover,  '-b');
    text(min(score_x_time_axis), 1.4*av_total_turnover,        num2str(av_total_turnover));
    title('Segment statistics');
    xlabel(xLabelText);
    ylabel(yScoreLabelText);
    legend('Turnover lamellipodium', 'Turnover lamella','Location','Best');
    
    
    hgsave(h_seg_av1,[figures_dir 'segment_av_turnover.fig']);
    print(h_seg_av1, [figures_dir 'segment_av_turnover.eps'],'-depsc2','-tiff');
    print(h_seg_av1, [figures_dir 'segment_av_turnover.tif'],'-dtiff');   
    
    h_seg_av2 = figure;
    subplot(2,1,1);
    plot(score_x_time_axis,  av_poly_time , 'b');
    hold on
    plot(score_x_time_axis,  av_shift_poly_time , '--b');
    plot(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_poly,  '-b');
    text(min(score_x_time_axis), 1.4*av_total_poly,      num2str(av_total_poly));
    
    ylabel(yScoreLabelTextplus);
    legend('Poly lamellipodium', 'Poly lamella', 'Location', 'Best');
    
    title('Segment statistics');
    
    subplot(2,1,2);
    plot(score_x_time_axis,  av_depoly_time , 'b');
    hold on
    plot(score_x_time_axis,  av_shift_depoly_time , '--b');
    plot(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_depoly,  '-b');
    text(min(score_x_time_axis), 1.4*av_total_depoly,      num2str(av_total_depoly));
    xlabel(xLabelText);
    ylabel(yScoreLabelTextminus);
    legend('De-poly lamellipodium', 'De-poly lamella','Location','Best');
    

    
    hgsave(h_seg_av2,[figures_dir 'segment_av_poly-depoly.fig']);
    print(h_seg_av2, [figures_dir 'segment_av_poly-depoly.eps'],'-depsc2','-tiff');
    print(h_seg_av2, [figures_dir 'segment_av_poly-depoly.tif'],'-dtiff');   
end




%%%% Plot only the protrusion and retrograde flow
if DO_PROT && ~DO_FLOW
    h_red_seg_av = figure;
    plot(prot_x_time_axis, av_protrusion_n_time,'g');
    hold on
    legend('Protrusion', 'Location','Best');
    plot(prot_x_time_axis, ones(size(prot_x_time_axis)).*av_total_prot_n,'-g');
    text(min(prot_x_time_axis), 1.4*av_total_prot_n, num2str(av_total_prot_n));
    title('Segment statistics');
    xlabel(xLabelText);
    ylabel(yLabelText);
elseif DO_PROT && DO_FLOW
    h_red_seg_av = figure;
    plot(prot_x_time_axis, av_protrusion_n_time,'g');
    hold on
    plot(flow_x_time_axis, av_retrograde_flow_n_time,       'r');
    plot(flow_x_time_axis, av_shift_retrograde_flow_n_time,  ':r');
    legend('Protrusion','Retrograde flow','Shifted retrograde flow', 'Location','Best');
    plot(prot_x_time_axis, ones(size(prot_x_time_axis)).*av_total_prot_n,'-g');
    text(min(prot_x_time_axis), 1.4*av_total_prot_n,             num2str(av_total_prot_n));
    plot(flow_x_time_axis, ones(size(flow_x_time_axis)).*av_total_retro_flow_n,'-r');
    plot(flow_x_time_axis, ones(size(flow_x_time_axis)).*av_total_shift_retro_flow_n,':r');
    text(min(flow_x_time_axis), 1.4*av_total_retro_flow_n,       num2str(av_total_retro_flow_n));
    text(min(flow_x_time_axis), 1.4*av_total_shift_retro_flow_n, num2str(av_total_shift_retro_flow_n));
    title(['Segment statistics']);
    xlabel(xLabelText);
    ylabel(yLabelText);
end

if DO_SCORES && DO_PROT
    % Set the second y axis for scores
%     ax1 = gca;
%     set(ax1,'box','off');
%     xlimits1 = get(ax1,'XLim');
%     ylimits1 = get(ax1,'YLim');
%     ylimits2(1) =     min(av_scores_time);
%     ylimits2(2) =     max(av_scores_time);
%     if abs(ylimits1(1) / ylimits1(2)) < abs(ylimits2(1) / ylimits2(2))
%         ylimits2(2) = (ylimits1(2) / ylimits1(1)) * ylimits2(1);
%     else
%         ylimits2(1) = (ylimits1(1) / ylimits1(2)) * ylimits2(2);
%     end
%     ax2 = axes('position',get(ax1,'position'));
%     set(ax2,'XAxisLocation','bottom','YAxisLocation','right','color','no');
%     axis(ax2,[xlimits1(1) xlimits1(2) ylimits2(1)  ylimits2(2)]);

    % plot average scores
    
    
    
%     line(score_x_time_axis ,av_scores_time ,'Color','b');
%     line(score_x_time_axis ,av_shift_scores_time ,'Color','b','LineStyle','--');
%     line(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_turnover, 'Color', 'b');
%     text(min(score_x_time_axis), 1.4*av_total_turnover,        num2str(av_total_turnover));
%     ylabel(yScoreLabelText);
    
    figure
    line(score_x_time_axis ,av_turnover_time ,'Color','b');
    line(score_x_time_axis ,av_shift_turnover_time ,'Color','b','LineStyle','--');
    line(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_turnover, 'Color', 'b');
    text(min(score_x_time_axis), 1.4*av_total_turnover,        num2str(av_total_turnover));
    title(['Segment statistics']);
    xlabel(xLabelText);
    ylabel(yScoreLabelText);
    legend('Turnover in Lamellipodium','Turnover in Lamella','Location','Best');
    
    figure
    line(score_x_time_axis ,av_poly_time ,'Color','b');
    line(score_x_time_axis ,av_shift_poly_time ,'Color','b','LineStyle','--');
    line(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_poly, 'Color', 'b');
    text(min(score_x_time_axis), 1.4*av_total_poly,        num2str(av_total_poly));
    title(['Segment statistics']);
    xlabel(xLabelText);
    ylabel(yScoreLabelTextplus);
    legend('Poly in Lamellipodium','Poly in Lamella','Location','Best');
    
    figure
    line(score_x_time_axis ,av_depoly_time ,'Color','b');
    line(score_x_time_axis ,av_shift_depoly_time ,'Color','b','LineStyle','--');
    line(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_depoly, 'Color', 'b');
    text(min(score_x_time_axis), 1.4*av_total_depoly,        num2str(av_total_depoly));
    title(['Segment statistics']);
    xlabel(xLabelText);
    ylabel(yScoreLabelTextminus);
    legend('De-Poly in Lamellipodium','De-Poly in Lamella','Location','Best');
end

if DO_PROT
    hgsave(h_red_seg_av,[figures_dir 'segment_av_prot_retro.fig']); 
    print(h_red_seg_av, [figures_dir 'segment_av_prot_retro.eps'],'-depsc2','-tiff'); 
    print(h_red_seg_av, [figures_dir 'segment_av_prot_retro.tif'],'-dtiff');
end

if DO_SCORES && ~DO_PROT
    % plot average scores
    figure
    line(score_x_time_axis ,av_turnover_time ,'Color','b');
    line(score_x_time_axis ,av_shift_turnover_time ,'Color','b','LineStyle','--');
    line(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_turnover, 'Color', 'b');
    text(min(score_x_time_axis), 1.4*av_total_turnover,        num2str(av_total_turnover));
    title(['Segment statistics']);
    xlabel(xLabelText);
    ylabel(yScoreLabelText); 
    legend('Turnover in Lamellipodium','Turnover in Lamella','Location','Best');
    
    figure
    line(score_x_time_axis ,av_poly_time ,'Color','b');
    line(score_x_time_axis ,av_shift_poly_time ,'Color','b','LineStyle','--');
    line(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_poly, 'Color', 'b');
    text(min(score_x_time_axis), 1.4*av_total_poly,        num2str(av_total_poly));
    title('Segment statistics');
    xlabel(xLabelText);
    ylabel(yScoreLabelTextplus);
    legend('Poly in Lamellipodium','Poly in Lamella','Location','Best');
    
    figure
    line(score_x_time_axis ,av_depoly_time ,'Color','b');
    line(score_x_time_axis ,av_shift_depoly_time ,'Color','b','LineStyle','--');
    line(score_x_time_axis,  ones(size(score_x_time_axis)).*av_total_depoly, 'Color', 'b');
    text(min(score_x_time_axis), 1.4*av_total_depoly,        num2str(av_total_depoly));
    title('Segment statistics');
    xlabel(xLabelText);
    ylabel(yScoreLabelTextminus);
    legend('De-Poly in Lamellipodium','De-Poly in Lamella','Location','Best');
    
end

if DO_ACTIVITY_1
    h_activity_1_seg_av = figure;
    plot(act_x_time_axis{1}, av_activity_time(1,:),'m');
    hold on
    if WINDOWS_TYPE == 1
        plot(act_x_time_axis{1}, av_shift_activity_time(1,:),'m--');
    end
    legend('Activity 1', 'Shift Activity 1', 'Location','Best');
    
    % plot average activity
    plot(act_x_time_axis{1},  ones(size(act_x_time_axis{1})).*av_total_activity(1), 'm');
    text(min(act_x_time_axis{1}), av_total_activity(1),        num2str(av_total_activity(1)));
    if WINDOWS_TYPE == 1
    plot(act_x_time_axis{1},  ones(size(act_x_time_axis{1})).*av_total_shift_activity(1),  'm--');
    text(min(act_x_time_axis{1}), av_total_shift_activity(1),        num2str(av_total_shift_activity(1)));
    end
    xlim([act_x_time_axis{1}(1) act_x_time_axis{1}(end)]);
    
    title(['Segment statistics']);   
    xlabel(xLabelText);
    ylabel('Intensity per pixel');     

    hgsave(h_activity_1_seg_av,[figures_dir 'segment_av_activity_1.fig']); 
    print(h_activity_1_seg_av, [figures_dir 'segment_av_activity_1.eps'],'-depsc2','-tiff'); 
    print(h_activity_1_seg_av, [figures_dir 'segment_av_activity_1.tif'],'-dtiff');  
end
if DO_ACTIVITY_2
    h_activity_2_seg_av = figure;
    plot(act_x_time_axis{2}, av_activity_time(2,:),'m');
    hold on
    if WINDOWS_TYPE == 1
        plot(act_x_time_axis{2}, av_shift_activity_time(2,:),'m--');
    end
    legend('Activity 2', 'Shift Activity 2', 'Location','Best');
    
    % plot average activity
    plot(act_x_time_axis{2},  ones(size(act_x_time_axis{2})).*av_total_activity(2), 'm');
    text(min(act_x_time_axis{2}), av_total_activity(2),        num2str(av_total_activity(2)));
    if WINDOWS_TYPE == 1
    plot(act_x_time_axis{2},  ones(size(act_x_time_axis{2})).*av_total_shift_activity(2),  'm--');
    text(min(act_x_time_axis{2}), av_total_shift_activity(2),        num2str(av_total_shift_activity(2)));
    end
    xlim([act_x_time_axis{2}(1) act_x_time_axis{2}(end)]);
    
    title('Segment statistics');
    xlabel(xLabelText);    
    ylabel('Intensity per pixel');   
    
    hgsave(h_activity_2_seg_av,[figures_dir 'segment_av_activity_2.fig']); 
    print(h_activity_2_seg_av, [figures_dir 'segment_av_activity_2.eps'],'-depsc2','-tiff'); 
    print(h_activity_2_seg_av, [figures_dir 'segment_av_activity_2.tif'],'-dtiff');  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot the angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DO_PROT
    h_angle = figure;
end
if DO_PROT && DO_FLOW
    plot(flow_x_time_axis, av_retrograde_flow_angle_time, '-sr')
    hold on
    plot(flow_x_time_axis, av_protrusion_angle_time, '-db')
    xlabel('Image #');
    ylabel('Angle [degree]');
    title('Segment statistics');
    legend('Retrograde flow - normal direction', 'Protrusion - normal direction','Location','Best');
elseif DO_PROT
    plot(prot_x_time_axis, av_protrusion_angle_time, '-db')
    xlabel('Image #');
    ylabel('Angle [degree]');
    title('Segment statistics');
    legend('Protrusion - normal direction', 'Location','Best');
end
if  DO_PROT
    hgsave(h_angle,[figures_dir 'angle_time_av.fig']); 
    print(h_angle, [figures_dir 'angle_time_av.eps'],'-depsc2','-tiff'); 
    print(h_angle, [figures_dir 'angle_time_av.tif'],'-dtiff');  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    Plot comparison of normal direction values %%%%%%%%%%
%%%%%%%%%    and non normal direction values            %%%%%%%%%%
x_axis = START_SEG:SEG_NR - END_SEG;
if DO_PROT
    figure
end
if DO_PROT && DO_FLOW
    %x_axis = START_SEG:SEG_NR - END_SEG;
    plot(x_axis, av_protrusion_n_seg,'-ob');
    hold on
    plot(x_axis, av_protrusion_seg,'-xb');

    plot(x_axis, av_retrograde_flow_n_seg,'-or');
    plot(x_axis, av_retrograde_flow_seg,'-xr');

    title('Average (normal) Protrusion in each segment');
    legend('normal protrusion', 'protrusion', 'normal retrograde flow', 'retrograde flow', 'Location', 'Best');
    line([START_SEG length(av_protrusion_seg)], [av_total_prot_n av_total_prot_n],'Color','b');
    line([START_SEG length(av_protrusion_seg)], [av_total_prot av_total_prot],'Color','b');   
    line([START_SEG length(av_protrusion_seg)], [av_total_retro_flow_n  av_total_retro_flow_n ],'Color','r');
    line([START_SEG length(av_protrusion_seg)], [av_total_retro_flow av_total_retro_flow],'Color','r');    
    text(START_SEG+1, 1.2*av_total_prot_n,        num2str(av_total_prot_n));
    text(START_SEG+1, 1.2*av_total_prot,          num2str(av_total_prot));
    text(START_SEG+1, 1.2*av_total_retro_flow_n,  num2str(av_total_retro_flow_n));
    text(START_SEG+1, 1.2*av_total_retro_flow,    num2str(av_total_retro_flow)); 
    
    xlabel('Segment #');
    ylabel('Velocity [pixel]');
elseif DO_PROT
    %x_axis = START_SEG:SEG_NR - END_SEG;
    plot(x_axis, av_protrusion_n_seg,'-ob');
    hold on
    plot(x_axis, av_protrusion_seg,'-xb');
    line([START_SEG length(av_protrusion_seg)], [av_total_prot_n av_total_prot_n],'Color','b');
    line([START_SEG length(av_protrusion_seg)], [av_total_prot av_total_prot],'Color','b');
    text(START_SEG+1, 1.7*av_total_prot_n,        num2str(av_total_prot_n));
    text(START_SEG+1, 1.04*av_total_prot,          num2str(av_total_prot));
    title('Average (normal) Protrusion in time');
    legend('normal protrusion', 'protrusion', 'Location', 'Best');
    xlabel('Segment #');
    ylabel('Velocity [pixel]');
end

if DO_PROT
    
    h_seg_av_disp = figure;
    plot(x_axis, av_displacement_n_seg,'-ob');
    hold on
    plot(x_axis, av_displacement_seg,'-xb');
    line([START_SEG length(av_displacement_seg)], [av_total_disp_n av_total_disp_n],'Color','b');
    line([START_SEG length(av_displacement_seg)], [av_total_disp av_total_disp],'Color','b');
    text(START_SEG+1, 1.7*av_total_disp_n,        num2str(av_total_disp_n));
    text(START_SEG+1, 1.04*av_total_disp,          num2str(av_total_disp));
    title('Average (normal) displacement in each segment');
    legend('normal displacement', 'displacement', 'Location', 'Best');
    xlabel('Segment #');
    ylabel(yLabelTextDisplacement);
    hgsave(h_seg_av_disp,[figures_dir 'seg_av_displacement.fig']);
    print(h_seg_av_disp, [figures_dir 'seg_av_displacement.eps'],'-depsc2','-tiff');
    print(h_seg_av_disp, [figures_dir 'seg_av_displacement.tif'],'-dtiff');    
    
    h_time_av_disp = figure;
    plot(prot_x_time_axis, av_displacement_n_time,'-ob');
    
    hold on
    plot(prot_x_time_axis, av_displacement_time,'-xb');
    line([FIRST_TIME length(av_displacement_time)*TIME_INTERVAL], [av_total_disp_n av_total_disp_n],'Color','b');
    line([FIRST_TIME length(av_displacement_time)*TIME_INTERVAL], [av_total_disp av_total_disp],'Color','b');
    text(START_SEG+1, 1.7*av_total_disp_n,        num2str(av_total_disp_n));
    text(START_SEG+1, 1.04*av_total_disp,          num2str(av_total_disp));
    title('Average (normal) displacement in time');
    legend('normal displacement', 'displacement', 'Location','Best');
    xlabel(xLabelText);
    ylabel(yLabelTextDisplacement);
    hgsave(h_time_av_disp,[figures_dir 'time_av_displacement.fig']);
    print(h_time_av_disp, [figures_dir 'time_av_displacement.eps'],'-depsc2','-tiff');
    print(h_time_av_disp, [figures_dir 'time_av_displacement.tif'],'-dtiff');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Plot the average number of samples  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%  per segment and time step    %%%%%%%%%%%%%%%%%%%
if DO_SCORES || DO_FLOW
    h_time_av_num_seg = figure;
    if DO_SCORES
        plot(x_axis, av_turnover_sample_num_seg ,'-xb');
        hold on
        plot(x_axis, av_shift_turnover_sample_num_seg ,'-sb');
        plot(x_axis, av_poly_sample_num_seg ,'->g');
        plot(x_axis, av_shift_poly_sample_num_seg ,'-sg');
        plot(x_axis, av_depoly_sample_num_seg ,'-<r');
        plot(x_axis, av_shift_depoly_sample_num_seg ,'-sr');
       
    end
    hold on
    if DO_FLOW
        plot(x_axis, av_num_seg_velocities_seg,'-xk');
        if WINDOWS_TYPE == 1
        plot(x_axis, av_num_seg_shift_velocities_seg, '-sk');
        end
    end
    if DO_SCORES && DO_FLOW
        legend('no. turnover', 'no. shift turnover', 'no. poly',...
            'no. shift poly','no. de-poly',...
            'no. shift de-poly','no. velocity vectors', 'no. shift velocity vectors','Location','Best');
    elseif DO_SCORES
        legend('no. turnover', 'no. shift turnover','no. poly',...
            'no. shift poly','no. de-poly','no. shift de-poly','Location','Best');
    else
        legend('no. velocity vectors', 'no. shift velocity vectors','Location','Best');
    end
    title('Ave. number of poly/depoly/turnover events per segment');
    %legend('# scores', '# shift scores', '# velocity vectors', '# shift velocity vectors');
    xlabel('Segment number');
    ylabel('counts');
    clear x_axis;
    hgsave(h_time_av_num_seg,[av_data_dir 'time_av_num_poly-depoly-turnover.fig']);
    print(h_time_av_num_seg, [av_data_dir 'time_av_num_poly-depoly-turnover.eps'],'-depsc2','-tiff');
    print(h_time_av_num_seg, [av_data_dir 'time_av_num_poly-depoly-turnover.tif'],'-dtiff');
    
    h_seg_av_num_time_score = figure;
    if DO_SCORES
        av_turnover_sample_num_time_plot = av_turnover_sample_num_time(:,2:end);
        plot(score_x_time_axis, av_turnover_sample_num_time_plot ,'-xb');
        hold on
        av_shift_turnover_sample_num_time_plot = av_shift_turnover_sample_num_time(:,2:end);
        plot(score_x_time_axis, av_shift_turnover_sample_num_time_plot ,'-sb');
        av_poly_sample_num_time_plot = av_poly_sample_num_time(:,2:end);
        plot(score_x_time_axis, av_poly_sample_num_time_plot ,'->g');
        av_shift_poly_sample_num_time_plot = av_shift_poly_sample_num_time(:,2:end);
        plot(score_x_time_axis, av_shift_poly_sample_num_time_plot ,'-sg');
        av_depoly_sample_num_time_plot = av_depoly_sample_num_time(:,2:end);
        plot(score_x_time_axis, av_depoly_sample_num_time_plot ,'-<r');
        av_shift_depoly_sample_num_time_plot = av_shift_depoly_sample_num_time(:,2:end);
        plot(score_x_time_axis, av_shift_depoly_sample_num_time_plot ,'-sr');
       
        legend('no. turnover', 'no. shift turnover', 'no. poly',...
            'no. shift poly','no. de-poly',...
            'no. shift de-poly','Location','Best');
        title('Ave. number of poly/depoly/turnover events and velocities at each time step');
        xlabel(xLabelText);
        ylabel('counts');
    hgsave(h_seg_av_num_time_score,[av_data_dir 'seg_av_num_poly-depoly-turnover.fig']);
    print(h_seg_av_num_time_score, [av_data_dir 'seg_av_num_poly-depoly-turnover.eps'],'-depsc2','-tiff');
    print(h_seg_av_num_time_score, [av_data_dir 'seg_av_num_poly-depoly-turnover.tif'],'-dtiff');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Plot the total number of turnover, poly, depoly %%%%%%%%%
%%%%%%%%%%%%%%%%  events in each segment and each time point  %%%%%%%%%%%%%
if DO_SCORES
    x_axis = START_SEG:SEG_NR - END_SEG;
%     h_time_total_num_seg = figure;
%     
%         plot(x_axis, total_turnover_sample_num_seg ,'-xb');
%         hold on
%         plot(x_axis, total_shift_turnover_sample_num_seg ,'-sb');
%         plot(x_axis, total_poly_sample_num_seg ,'->g');
%         plot(x_axis, total_shift_poly_sample_num_seg ,'-sg');
%         plot(x_axis, total_depoly_sample_num_seg ,'-<r');
%         plot(x_axis, total_shift_depoly_sample_num_seg ,'-sr');
% 
%         legend('tot. no. turnover', 'tot. no. shift turnover', 'tot. no. poly',...
%             'tot. no. shift poly','tot. no. de-poly',...
%             'tot. no. shift de-poly','Location','Best');
%     
%     title('Total number of samples per segment');
%     %legend('# scores', '# shift scores', '# velocity vectors', '# shift velocity vectors');
%     xlabel('Segment number');
%     ylabel('counts');
%     hgsave(h_time_total_num_seg,[av_data_dir 'time_total_num_poly-depoly-turnover.fig']);
%     print(h_time_total_num_seg, [av_data_dir 'time_total_num_poly-depoly-turnover.eps'],'-depsc2','-tiff');
%     print(h_time_total_num_seg, [av_data_dir 'time_total_num_poly-depoly-turnover.tif'],'-dtiff');
    
    h_seg_total_num_time_score = figure;
    
        total_turnover_sample_num_time_plot = total_turnover_sample_num_time(:,2:end);
        plot(score_x_time_axis, total_turnover_sample_num_time_plot ,'-xb');
        hold on
        total_shift_turnover_sample_num_time_plot = total_shift_turnover_sample_num_time(:,2:end);
        plot(score_x_time_axis, total_shift_turnover_sample_num_time_plot ,'-sb');
        total_poly_sample_num_time_plot = total_poly_sample_num_time(:,2:end);
        plot(score_x_time_axis, total_poly_sample_num_time_plot ,'->g');
        total_shift_poly_sample_num_time_plot = total_shift_poly_sample_num_time(:,2:end);
        plot(score_x_time_axis, total_shift_poly_sample_num_time_plot ,'-sg');
        total_depoly_sample_num_time_plot = total_depoly_sample_num_time(:,2:end);
        plot(score_x_time_axis, total_depoly_sample_num_time_plot ,'-<r');
        total_shift_depoly_sample_num_time_plot = total_shift_depoly_sample_num_time(:,2:end);
        plot(score_x_time_axis, total_shift_depoly_sample_num_time_plot ,'-sr');
       
        legend('tot. no. turnover', 'tot. no. shift turnover', 'tot. no. poly',...
            'tot. no. shift poly','tot. no. de-poly',...
            'tot. no. shift de-poly','Location','Best');
        title('Tot. number of poly/depoly/turnover events at each time step');
        xlabel(xLabelText);
        ylabel('counts');
    hgsave(h_seg_total_num_time_score,[av_data_dir 'seg_tot_num_poly-depoly-turnover.fig']);
    print(h_seg_total_num_time_score, [av_data_dir 'seg_tot_num_poly-depoly-turnover.eps'],'-depsc2','-tiff');
    print(h_seg_total_num_time_score, [av_data_dir 'seg_tot_num_poly-depoly-turnover.tif'],'-dtiff');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     if DO_FLOW
%         h_seg_av_num_time_flow = figure;
%         plot(x_time_axis, av_num_seg_velocities_time,'-xk');
%         hold on;
%         if WINDOWS_TYPE == 1
%         plot(x_time_axis, av_num_seg_shift_velocities_time, '-sk');
%         end
%         legend('no. velocity vectors', 'no. shift velocity vectors','Location','Best');
%         title('Ave. number of velocity vectors in Lamella and Lamellipodium in time');
%         xlabel(xLabelText);
%         ylabel('counts');
%     hgsave(h_seg_av_num_time_flow,[av_data_dir 'seg_av_num_velocities.fig']);
%     print(h_seg_av_num_time_flow, [av_data_dir 'seg_av_num_velocities.eps'],'-depsc2','-tiff');
%     print(h_seg_av_num_time_flow, [av_data_dir 'seg_av_num_velocities.tif'],'-dtiff');
%     end
%     if DO_SCORES && DO_FLOW
%         legend('# turnover', '# shift turnover', '# poly',...
%             '# shift poly','# de-poly',...
%             '# shift de-poly','# velocity vectors', '# shift velocity vectors','Location','Best');
%     elseif DO_SCORES
%         legend('# turnover', '# shift turnover','# poly',...
%             '# shift poly','# de-poly','# shift de-poly','Location','Best');
%     else
%         legend('# velocity vectors', '# shift velocity vectors','Location','Best');
%     end
    
    %legend('# scores', '# shift scores', '# velocity vectors', '# shift velocity vectors');
%     xlabel(xLabelText);
%     ylabel('#');
%     clear x_axis;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if WINDOWS_TYPE == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Convert to images   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (TOTAL_TIME_STEPS - FIRST_TIME) > 5
    if DO_SCORES
        img_turnover = imresize(score,IMAGE_STRETCH, 'nearest');
        img_shift_turnover = imresize(shift_score,IMAGE_STRETCH, 'nearest');
        img_poly = imresize(score_plus,IMAGE_STRETCH, 'nearest');
        img_shift_poly = imresize(shift_score_plus,IMAGE_STRETCH, 'nearest');
        img_depoly = imresize(score_minus,IMAGE_STRETCH, 'nearest');
        img_shift_depoly = imresize(shift_score_minus,IMAGE_STRETCH, 'nearest');
    end
    if DO_PROT
        img_protrusion  = imresize(protrusion_normal,IMAGE_STRETCH, 'nearest');
    end
    if DO_FLOW
        img_retrograde_flow         = imresize(retrograde_flow_normal,IMAGE_STRETCH, 'nearest');
        img_shift_retrograde_flow   = imresize(shift_retrograde_flow_normal,IMAGE_STRETCH, 'nearest');
        img_av_network_velocity_direction_cor = imresize(av_network_velocity_direction_cor,IMAGE_STRETCH, 'nearest');
    end
    if DO_ACTIVITY_1
        img_activity_1 = imresize(squeeze(activity(1,:,:)),IMAGE_STRETCH, 'nearest');
        img_shift_activity_1 = imresize(squeeze(shift_activity(1,:,:)),IMAGE_STRETCH, 'nearest');
    end
    if DO_ACTIVITY_2
        img_activity_2 = imresize(squeeze(activity(2,:,:)),IMAGE_STRETCH, 'nearest');
        img_shift_activity_2 = imresize(squeeze(shift_activity(2,:,:)),IMAGE_STRETCH, 'nearest');
    end    
    %%%%%%%%%%%   Filter              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filter_kernel = fspecial('gaussian',[IMG_GAUSS_W IMG_GAUSS_W], IMG_GAUSS_SIG);
    if DO_SCORES
        [img_turnover_f, dum] = Gauss2DBorder(img_turnover, IMG_GAUSS_SIG);
        save([merg_dir 'img_turnover.mat'], 'img_turnover', 'img_turnover_f');
        [img_shift_turnover_f, dum] = Gauss2DBorder(img_shift_turnover, IMG_GAUSS_SIG);
        save([merg_dir 'img_shift_turnover.mat'], 'img_shift_turnover', 'img_shift_turnover_f'); 
        [img_poly_f, dum] = Gauss2DBorder(img_poly, IMG_GAUSS_SIG);
        save([merg_dir 'img_poly.mat'], 'img_poly', 'img_poly_f');
        [img_shift_poly_f, dum] = Gauss2DBorder(img_shift_poly, IMG_GAUSS_SIG);
        save([merg_dir 'img_shift_poly.mat'], 'img_shift_poly', 'img_shift_poly_f');
        [img_depoly_f, dum] = Gauss2DBorder(img_depoly, IMG_GAUSS_SIG);
        save([merg_dir 'img_depoly.mat'], 'img_depoly', 'img_depoly');
        [img_shift_depoly_f, dum] = Gauss2DBorder(img_shift_depoly, IMG_GAUSS_SIG);
        save([merg_dir 'img_shift_depoly.mat'], 'img_shift_depoly', 'img_shift_depoly_f');
    end
    if DO_PROT
        [img_protrusion_f, dum] = Gauss2DBorder(img_protrusion, IMG_GAUSS_SIG);
        save([merg_dir 'img_prot.mat'], 'img_protrusion',  'img_protrusion_f');
    end
    if DO_FLOW
        [img_retrograde_flow_f, dum] = Gauss2DBorder(img_retrograde_flow, IMG_GAUSS_SIG);
        [img_shift_retrograde_flow_f, dum] = Gauss2DBorder(img_shift_retrograde_flow, IMG_GAUSS_SIG);
        save([merg_dir 'img_retrograde.mat'], 'img_retrograde_flow',  'img_shift_retrograde_flow',...
            'img_retrograde_flow_f', 'img_shift_retrograde_flow_f');
        [img_av_network_velocity_direction_cor_f, dum] = Gauss2DBorder(img_av_network_velocity_direction_cor, IMG_GAUSS_SIG);
        save([merg_dir 'img_av_network_velocity_direction_cor_f.mat'],'img_av_network_velocity_direction_cor',...
            'img_av_network_velocity_direction_cor_f');
    end
    if DO_ACTIVITY_1
        [img_activity_1_f, dum] = Gauss2DBorder(img_activity_1, IMG_GAUSS_SIG);
        [img_shift_activity_1_f, dum] = Gauss2DBorder(img_shift_activity_1, IMG_GAUSS_SIG);
        save([merg_dir 'img_activity_1.mat'], 'img_activity_1',   'img_activity_1_f');
        save([merg_dir 'img_shift_activity_1.mat'], 'img_shift_activity_1',   'img_shift_activity_1_f');
    end
    if DO_ACTIVITY_2
        [img_activity_2_f, dum] = Gauss2DBorder(img_activity_2, IMG_GAUSS_SIG);
        [img_shift_activity_2_f, dum] = Gauss2DBorder(img_shift_activity_2, IMG_GAUSS_SIG);
        save([merg_dir 'img_activity_2.mat'], 'img_activity_2',   'img_activity_2_f');
        save([merg_dir 'img_shift_activity_2.mat'], 'img_shift_activity_2',   'img_shift_activity_2_f');
    end
    
%%%% modifies on 7.26.07 by Mohsen to make score, activity1,2 and 
%%%% protrusion figures saved in \merge\figures\ directory show time stamp
%%%% and proper labling.

        No_Ticks = 10;
        iptsetpref('ImshowAxesVisible','on');

        if DO_SCORES
        h_score_img = figure;
        imshow(img_turnover_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['Lamellipodium ' yScoreLabelText]);
        xlabel(xLabelText);
        ylabel('Segment #');
        c = redGreenColorMap;
        colormap(c);
        MAX_COLOR_VAL_SCORES  = abs(robust_av_turnover(1)+0.5*...
             robust_av_turnover(2)*score_x_time_axis(end)) + 2*robust_sigma_turnover;
        caxis([-MAX_COLOR_VAL_SCORES MAX_COLOR_VAL_SCORES]);
      
        colorbar;
        hgsave(h_score_img,[figures_dir 'Turnover_activity_img.fig']);
        print(h_score_img, [figures_dir 'Turnover_activity_img.eps'],'-depsc2','-tiff');
        print(h_score_img, [figures_dir 'Turnover_activity_img.tif'],'-dtiff');
        
        h_shift_score_img = figure;
        imshow(img_shift_turnover_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['lamella ' yScoreLabelText]);
        xlabel(xLabelText);
        ylabel('Segment #');
        c = redGreenColorMap;
        colormap(c);
        colorbar;
        MAX_COLOR_VAL_SCORES = abs(robust_av_shift_turnover(1)+...
            0.5*robust_av_shift_turnover(2)*score_x_time_axis(end)) +...
            2*robust_sigma_shift_turnover;
        caxis([-MAX_COLOR_VAL_SCORES MAX_COLOR_VAL_SCORES]);
        colorbar;
        hgsave(h_shift_score_img,[figures_dir 'shift_Turnover_activity_img.fig']);
        print(h_shift_score_img, [figures_dir 'shift_Turnover_activity_img.eps'],'-depsc2','-tiff');
        print(h_shift_score_img, [figures_dir 'shift_Turnover_activity_img.tif'],'-dtiff'); 
        
        h_score_img_plus = figure;
        imshow(img_poly_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['Lamellipodium ' yScoreLabelTextplus]);
        xlabel(xLabelText);
        ylabel('Segment #');
        c = jet;
        colormap(c);
        MAX_COLOR_VAL_SCORES_plus  = abs(robust_av_poly(1)+0.5*...
             robust_av_poly(2)*score_x_time_axis(end)) + 2*robust_sigma_poly;
        caxis([-MAX_COLOR_VAL_SCORES_plus MAX_COLOR_VAL_SCORES_plus]);
      
        colorbar;
        hgsave(h_score_img_plus,[figures_dir 'poly_activity_img.fig']);
        print(h_score_img_plus, [figures_dir 'poly_activity_img.eps'],'-depsc2','-tiff');
        print(h_score_img_plus, [figures_dir 'poly_activity_img.tif'],'-dtiff');
        
        h_shift_score_img_plus = figure;
        imshow(img_shift_poly_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['lamella ' yScoreLabelTextplus]);
        xlabel(xLabelText);
        ylabel('Segment #');
        c = jet;
        colormap(c);
        MAX_COLOR_VAL_SCORES = abs(robust_av_shift_poly(1)+... 
            0.5*robust_av_shift_poly(2)*score_x_time_axis(end)) +...
            2*robust_sigma_shift_poly;
        caxis([-MAX_COLOR_VAL_SCORES MAX_COLOR_VAL_SCORES]);
        colorbar;
        hgsave(h_shift_score_img_plus,[figures_dir 'shift_poly_activity_img.fig']);
        print(h_shift_score_img_plus, [figures_dir 'shift_poly_activity_img.eps'],'-depsc2','-tiff');
        print(h_shift_score_img_plus, [figures_dir 'shift_poly_activity_img.tif'],'-dtiff'); 
        
        h_score_img_minus = figure;
        imshow(img_depoly_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['Lamellipodium ' yScoreLabelTextminus]);
        xlabel(xLabelText);
        ylabel('Segment #');
        c = jet;
        colormap(c);
        MAX_COLOR_VAL_SCORES_minus  = abs(robust_av_depoly(1)+0.5*...
             robust_av_depoly(2)*score_x_time_axis(end)) + 2*robust_sigma_depoly;
        caxis([-MAX_COLOR_VAL_SCORES_minus MAX_COLOR_VAL_SCORES_minus]);
      
        colorbar;
        hgsave(h_score_img_minus,[figures_dir 'depoly_activity_img.fig']);
        print(h_score_img_minus, [figures_dir 'depoly_activity_img.eps'],'-depsc2','-tiff');
        print(h_score_img_minus, [figures_dir 'depoly_activity_img.tif'],'-dtiff');
        
        h_shift_score_img_minus = figure;
        imshow(img_shift_depoly_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['lamella ' yScoreLabelTextminus]);
        xlabel(XlabelText);
        ylabel('Segment #');
        c = jet;
        colormap(c);
        MAX_COLOR_VAL_SCORES = abs(robust_av_shift_depoly(1)+... 
            0.5*robust_av_shift_depoly(2)*score_x_time_axis(end)) +...
            2*robust_sigma_shift_depoly;
        caxis([-MAX_COLOR_VAL_SCORES MAX_COLOR_VAL_SCORES]);
        colorbar;
        hgsave(h_shift_score_img_minus,[figures_dir 'shift_depoly_activity_img.fig']);
        print(h_shift_score_img_minus, [figures_dir 'shift_depoly_activity_img.eps'],'-depsc2','-tiff');
        print(h_shift_score_img_minus, [figures_dir 'shift_depoly_activity_img.tif'],'-dtiff');
        
    end

    if DO_ACTIVITY_1
        h_activity_img = figure;
        imshow(img_activity_1_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title('Activity 1 (Intensity/Pixel)');
        xlabel(xLabelText);
        ylabel('Segment #');
        colormap(hot);
        MAX_COLOR_VAL_ACTIVITY = abs(robust_av_activity{1}(1)+0.5*robust_av_activity{1}(2)*act_x_time_axis{1}(end)) + 5*robust_sigma_activity(1);
        caxis([0 MAX_COLOR_VAL_ACTIVITY]);
        colorbar;
        hgsave(h_activity_img,[figures_dir 'activity_1_img.fig']);
        print(h_activity_img, [figures_dir 'activity_1_img.eps'],'-depsc2','-tiff');
        print(h_activity_img, [figures_dir 'activity_1_img.tif'],'-dtiff');
        
        h_shift_activity_img = figure;
        imshow(img_shift_activity_1_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title('Activity 1 in Lamellum (Intensity/Pixel)');
        xlabel(xLabelText);
        ylabel('Segment #');
        colormap(hot);
        MAX_COLOR_VAL_ACTIVITY = abs(robust_av_shift_activity{1}(1)+0.5*robust_av_shift_activity{1}(2)*act_x_time_axis{1}(end)) + 5*robust_sigma_shift_activity(1);
        caxis([0 MAX_COLOR_VAL_ACTIVITY]);
        colorbar;
        hgsave(h_activity_img,[figures_dir 'shift_activity_1_img.fig']);
        print(h_activity_img, [figures_dir 'shift_activity_1_img.eps'],'-depsc2','-tiff');
        print(h_activity_img, [figures_dir 'shift_activity_1_img.tif'],'-dtiff');     
    end
    
    if DO_ACTIVITY_2
        h_activity_img = figure;
        imshow(img_activity_2_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        MAX_COLOR_VAL_ACTIVITY = abs(robust_av_activity{2}(1)+0.5*robust_av_activity{2}(2)*act_x_time_axis{2}(end)) + 5*robust_sigma_activity(2);
        caxis([0 MAX_COLOR_VAL_ACTIVITY]);
        colormap(hot);
        colorbar;
        title('Activity 2 (Intensity/Pixel) ');
        xlabel(xLabelText);
        ylabel('Segment #');
        hgsave(h_activity_img,[figures_dir 'activity_2_img.fig']);
        print(h_activity_img, [figures_dir 'activity_2_img.eps'],'-depsc2','-tiff');
        print(h_activity_img, [figures_dir 'activity_2_img.tif'],'-dtiff');
        
        h_shift_activity_img = figure;
        imshow(img_shift_activity_2_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title('Activity 2 in Lamellum (Intensity/Pixel)');
        xlabel(xLabelText);
        ylabel('Segment #');
        colormap(hot);
        MAX_COLOR_VAL_ACTIVITY= abs(robust_av_shift_activity{2}(1)+0.5*robust_av_shift_activity{2}(2)*act_x_time_axis{2}(end)) + 5*robust_sigma_shift_activity(2);
        caxis([0 MAX_COLOR_VAL_ACTIVITY]);
        colorbar;
        hgsave(h_activity_img,[figures_dir 'shift_activity_2_img.fig']);
        print(h_activity_img, [figures_dir 'shift_activity_2_img.eps'],'-depsc2','-tiff');
        print(h_activity_img, [figures_dir 'shift_activity_2_img.tif'],'-dtiff');     
    end
    
    if DO_PROT
        h_protrusion_img = figure;
        imshow(img_protrusion_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title('Activity 1 in Lamellum (Intensity/Pixel)');
        xlabel(xLabelText);
        ylabel('Segment #');
        title(['Protrusion' yLabelText]);
        colormap(jet);
        MAX_COLOR_VAL_PROT = abs(robust_av_protrusion(1)+0.5*robust_av_protrusion(2)*prot_x_time_axis(end)) + 3*robust_sigma_protrusion;
        caxis([-MAX_COLOR_VAL_PROT MAX_COLOR_VAL_PROT]);
        colorbar;
        hgsave(h_protrusion_img,[figures_dir 'prot_activity_img.fig']);
        print(h_protrusion_img, [figures_dir 'prot_activity_img.eps'],'-depsc2','-tiff');
        print(h_protrusion_img, [figures_dir 'prot_activity_img.tif'],'-dtiff');
    end

    if DO_FLOW
        
        h_retrograde_flow_img = figure;
        imshow(img_retrograde_flow_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['Retrograde flow in Lamellipodium' yLabelText]);
        xlabel(xLabelText);
        ylabel('Segment #');
        colormap(jet);
        MAX_COLOR_VAL_RETRO = abs(robust_av_retrograde_flow(1) + 0.5*robust_av_retrograde_flow(2)*flow_x_time_axis(end)) + 2*robust_sigma_retrograde_flow;
        caxis([-MAX_COLOR_VAL_RETRO MAX_COLOR_VAL_RETRO]);
        colorbar;
        hgsave(h_retrograde_flow_img,[figures_dir 'flow_activity_img.fig']);
        print(h_retrograde_flow_img, [figures_dir 'flow_activity_img.eps'],'-depsc2','-tiff');
        print(h_retrograde_flow_img, [figures_dir 'flow_activity_img.tif'],'-dtiff');
        
        h_shift_retrograde_flow_img = figure;
        imshow(img_shift_retrograde_flow_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['Retrograde flow in lamellum ' yLabelText]);
        xlabel(xLabelText);
        ylabel('Segment #');
        colormap(jet);
        caxis([-MAX_COLOR_VAL_RETRO MAX_COLOR_VAL_RETRO]);
        colorbar;
        hgsave(h_shift_retrograde_flow_img,[figures_dir 'shift_flow_activity_img.fig']);
        print(h_shift_retrograde_flow_img, [figures_dir 'shift_flow_activity_img.eps'],'-depsc2','-tiff');
        print(h_shift_retrograde_flow_img, [figures_dir 'shift_flow_activity_img.tif'],'-dtiff');

        
        
        h_network_velocity_direction_cor = figure;
        imshow(img_av_network_velocity_direction_cor_f);
        V = axis;
        XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
        XTickV = V(1,1):XTickV_STEP:V(1,2);
        YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
        YTickV = V(1,3):YTickV_STEP:V(1,4);
        set(gca,'XTick', XTickV);
        set(gca,'YTick', YTickV);
        XTickL = {};
        YTickL = {};
        for jj = 1:length(XTickV)
            XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
        end
        for jj = 1:length(YTickV)
            YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
        end
        set(gca,'XTickLabel',XTickL);
        set(gca,'YTickLabel',YTickL);
        title(['Retrograde flow direction correlation (scalar product)' ]);
        xlabel('Time');
        ylabel('Segment #');
        colormap(jet);
        caxis([-MAX_COLOR_VAL_RETRO MAX_COLOR_VAL_RETRO]);
        colorbar;
        hgsave(h_network_velocity_direction_cor,[figures_dir 'network_velocity_direction_cor_img.fig']);
        print(h_network_velocity_direction_cor, [figures_dir 'network_velocity_direction_cor_img.eps'],'-depsc2','-tiff');
        print(h_network_velocity_direction_cor, [figures_dir 'network_velocity_direction_cor_img.tif'],'-dtiff');        
    end
end
end  % if WINDOWS_TYPE == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose('all');










    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%  Create smart segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
        % create distance map
        dist_map = bwdist(~mask_img);
        % get the levels
        figure,imshow(mask_img);
        co = cell(50,1)
        hold on
        
        co_0 = contour(dist_map,[0.01 0.01])';
        [edge_sp_x(1), edge_sp_y(1)] = imPixelChainSpline(co_0, 'tolerance',15);  
        
        de = [15, 30, 45, 60];
        
        for i_level = 1:4
            co{i_level+1} = contour(dist_map,[de(i_level) de(i_level)])';
            co{i_level+1}(1:3,:)=[];
            co{i_level+1}(end-3:end,:)=[];
            hold on
            plot(co{i_level+1}(:,1),co{i_level+1}(:,2),'-');
            [edge_sp_x(i_level+1), edge_sp_y(i_level+1)] = imPixelChainSpline(co{i_level+1}, 'tolerance',15);
        end

        max_spline_par = edge_sp_x(1).knots(end);
        s_p = (max_spline_par-1)/SEG_NR;
        m_pos=1:s_p:max_spline_par;
        
        seg_pos(:,1) = m_pos';
        
        for i_level = 1:4
            [temp1, temp2, seg_pos(:,i_level+1)]=prGetDispNearest(edge_sp_x(i_level), edge_sp_y(i_level),...
                edge_sp_x(i_level+1), edge_sp_y(i_level+1), seg_pos(:,i_level)',...
                'tol', 20, 'robust_min', 20);

        end
        for i_level = 1:5
            seg_x(:,i_level) =  fnval(edge_sp_x(i_level),seg_pos(:,i_level));
            seg_y(:,i_level) =  fnval(edge_sp_y(i_level),seg_pos(:,i_level));
            plot(seg_x(:,i_level),seg_y(:,i_level),'rx');
        end
        a=1;
%
%
%         max_spline_par = edge_sp_x(1).knots(end);
%         s_p = (max_spline_par-1)/SEG_NR;
%         m_pos=1:s_p:max_spline_par;
%         
%         seg_pos(:,1) = m_pos;
%         for i_level = 1:50-1
%             if 0
%                 [temp1, temp2, i_pos]=prGetDispNearest(edge_sp_x(i_level), edge_sp_y(i_level),...
%                     edge_sp_x(i_level+1), edge_sp_y(i_level+1), m_pos,...
%                     'tol', 20, 'robust_min', 40);
%             elseif 0
%                 [temp1, temp2, i_pos]=prGetDispNormal(edge_sp_x(i_level), edge_sp_y(i_level),...
%                                                    edge_sp_x(i_level+1), edge_sp_y(i_level+1), m_pos);
%             else
%                                                    
%             end
%                                                    
%             m_pos = i_pos;
%             
%             seg_pos(:,i_level+1) = i_pos;
%         end
%         
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


