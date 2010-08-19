function prAlpha(varargin)

This function is replaced!!!!!!!
%
%
% SYNOPSIS      prAlpha
%
% INPUT                :         
% 
% OUTPUT               : 
%                           
% DEPENDENCES       prAlpha uses {    netAssemblyMaps
%                                     velocityMapsOrg                             
%                                       }
%
%                   prAlpha is used by { 
%                                           }
%
% Matthias Machacek 12/04/03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);
for i=1:l
    if strcmp(varargin(i),'pixel')
        PIXEL=varargin{i+1};  
    elseif strcmp(varargin(i),'bit_depth')
        BIT_DEPTH=varargin{i+1};          
    elseif strcmp(varargin(i),'prot_depth')
        PROT_DEPTH=varargin{i+1};  
    elseif strcmp(varargin(i),'time_interval')
        TIME_INTERVAL=varargin{i+1};               
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('PIXEL','var')
    %pixel size in nm
    PIXEL=67;
end
if ~exist('BIT_DEPTH','var')
    %the bit depth of the images
    BIT_DEPTH=16383;
end
if ~exist('PROT_DEPTH','var')
    %protrusion band mask width
    PROT_DEPTH=20;
end
if ~exist('TIME_INTERVAL','var')
    %time interval between the image frames [s]
    TIME_INTERVAL=10;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contr=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1 
    %%%%%%%%Project cell1
    % number of images  : 1-30 (30)
    % type of cell      :
    % drug              :
    dir_w =     'L:\projects\rho_protrusion\cell1\';
    img_dir =   'L:\projects\rho_protrusion\cell1\images\';
    %prot_dir =  'protrusion\';
    prot_dir =  'protrusion\';  
    img_name =  'plane01';
elseif 0   
    %%%%%%%%Project cell2
    % number of images  : 1-30 (30)
    % type of cell      :
    % drug              :
    dir_w=      'L:\projects\rho_protrusion\cell2\';
    img_dir =   'L:\projects\rho_protrusion\cell2\images\';    
    prot_dir =  'protrusion\';
    img_name =  'plane01';
end    


time_index          = 0;
FIRST_TIME          = 1;
TOTAL_TIME_STEPS    = 29;
PROT_DEPTH          = 10;
%this denotes the first segment used for building the average values
START_SEG   =1;
%this denotes how many less segments are used for the building the 
%average values
END_SEG     =2;
   
%offset [pixel] of the shifted observation region
%for the retrograde flow
SEG_SHIFT           = 40;

%color map maximal values for the activity maps
AUTO_COL_SCALE      =0;
MAX_COLOR_VAL_PROT  = 14;


%parameters for the smoothing of the activity images
IMG_GAUSS_W         = 20; %20
IMG_GAUSS_SIG       = 15;  %6
IMAGE_STRETCH       = 15;


PROCRUSTE           = 0;
WAVE_FIT            = 0;
SURF_ACTIVITY_PLOT  = 0;
ALL_DATA_CALC       = 0;
BOOTSTRAP           = 0;
IMG_DIFF            = 0;
VAR_SPACE_PLOT      = 0; %1
PROT_PLOT           = 1;


%results directory
dir_alpha = [dir_w 'pr_results' filesep];

%mask file names
firstfilename_mask=([dir_w prot_dir 'mask_' img_name '.tif']);
[filelist_mask]=getFileStackNames(firstfilename_mask);

% get the filelist filenames of the segment boundary masks
firstfilename_seg=([dir_w prot_dir  filesep 's_mask_' img_name '.dat']);
[filelist_seg]=getFileStackNames(firstfilename_seg);

% get the filelist filenames of the images
firstfilename_img=([img_dir img_name '.tif']);
[filelist_img]=getFileStackNames(firstfilename_img);

%file name containing the edge pixels
file_pixel_edge=[dir_w  prot_dir 'pixel_edge.dat'];

%file name containing the spline pixels
file_spline=[dir_w  prot_dir 'edge_spline.mat'];

%file name containing the averaged normal vectors
file_av_normal=[dir_w  prot_dir 'av_normal.dat'];

% get the filelist filenames containing the averaged protrusion vectors
file_av_prot=[dir_w  prot_dir 'av_prot.dat'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Write all the parameters to a file  %%%%%%%%%%%%%%%%%%
fid = fopen([dir_alpha 'parameters.dat'],'w+');
if fid == -1
    error('Could not create file in results directory. Be sure to put a "protrusion" directory into results directory!!');
end
fprintf(fid, 'Image file name = %s\n',          img_name);
fprintf(fid, 'Protrusion directory = %s\n',     prot_dir);
fprintf(fid, 'FIRST_TIME = %d\n',               FIRST_TIME);
fprintf(fid, 'TOTAL_TIME_STEPS = %d\n',         TOTAL_TIME_STEPS);
fprintf(fid, 'START_SEG = %d\n',                START_SEG);
fprintf(fid, 'END_SEG = %d\n',                  END_SEG);
fprintf(fid, 'MAX_COLOR_VAL_PROT = %d\n',       MAX_COLOR_VAL_PROT);
fprintf(fid, 'IMG_GAUSS_W = %d\n',              IMG_GAUSS_W);
fprintf(fid, 'IMG_GAUSS_SIG = %d\n',            IMG_GAUSS_SIG);
fprintf(fid, 'SEG_SHIFT = %d\n',                SEG_SHIFT);
fprintf(fid, 'PROT_DEPTH= %d\n',                PROT_DEPTH);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% open files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_pixel_edge= fopen(file_pixel_edge,'r');
if fid_pixel_edge == -1
    error('Could not open file pixel edge');
end

fid_av_normal = fopen(file_av_normal,'r');
if fid_av_normal == -1
    error('Could not open file averaged normals');
end
 
fid_av_prot = fopen(file_av_prot,'r');
if fid_av_prot == -1
    error('Could not open segment boundary coordinates');
end      

load(file_spline);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_pix     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
x_p_edge  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix]);
y_p_edge  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%wind the data to the desired starting time step
for time = 1: FIRST_TIME - 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the pixel edge at next time step%%%%%%%%%%%%
    n_pix_next     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
    x_p_edge_next  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_next]);
    y_p_edge_next  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_next]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the averaged normal vectors %%%%%%%%%%%%%%%%
    n_el_norm        = fscanf(fid_av_normal,'%g  ', [1 1]);
    spline_p_norm    = fscanf(fid_av_normal,'%g  ', [1 n_el_norm]);
    x_av_pos_normal  = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]);
    y_av_pos_normal  = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]); 
    x_av_normal      = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]);
    y_av_normal      = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the averaged protrusion along the edge %%%%%
    n_el_av_prot    = fscanf(fid_av_prot,'%g  ', [1 1]);
    spline_p_av_prot= fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot]);
    x_av_pos_prot   = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]);
    y_av_pos_prot   = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]);
    x_av_prot       = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]);
    y_av_prot       = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the protrusion along the edge %%%%%%%%%%%%%%
    n_el_prot       = fscanf(fid_prot,'%g  ', [1 1]);
    spline_p_prot   = fscanf(fid_prot,'%g  ', [1 n_el_prot]);
    x_pos_prot      = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]);
    y_pos_prot      = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]);
    x_prot          = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]);
    y_prot          = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




%the first time step refers to number of the files in the folder% 
%and NOT to the number in the file name!!!!!!!!!!!!   %%%%%%%%%%%

for time = FIRST_TIME : 1: TOTAL_TIME_STEPS
    time_index = time_index+1;
    %display the time step
    time
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fileName_img=char(filelist_img(time));
    image=imreadnd2(fileName_img,0,BIT_DEPTH);
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
    %%%%%%%%%%%%%%%% read the pixel edge at next time step%%%%%%%%%%%%
    n_pix_next     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
    x_p_edge_next  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_next]);
    y_p_edge_next  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_next]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the averaged normal vectors %%%%%%%%%%%%%%%%
    n_el_norm        = fscanf(fid_av_normal,'%g  ', [1 1]);
    spline_p_norm    = fscanf(fid_av_normal,'%g  ', [1 n_el_norm]);
    x_av_pos_normal  = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]);
    y_av_pos_normal  = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]); 
    x_av_normal      = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]);
    y_av_normal      = fscanf(fid_av_normal,'%g  ', [1 n_el_norm-1]);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the averaged protrusion along the edge %%%%%
    n_el_av_prot    = fscanf(fid_av_prot,'%g  ', [1 1]);
    spline_p_av_prot= fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot]);
    x_av_pos_prot   = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]);
    y_av_pos_prot   = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]);
    x_av_prot       = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]);
    y_av_prot       = fscanf(fid_av_prot,'%g  ', [1 n_el_av_prot-1]); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% read the protrusion along the edge %%%%%%%%%%%%%%
%     n_el_prot       = fscanf(fid_prot,'%g  ', [1 1]);
%     spline_p_prot   = fscanf(fid_prot,'%g  ', [1 n_el_prot]);
%     x_pos_prot      = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]);
%     y_pos_prot      = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]);
%     x_prot          = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]);
%     y_prot          = fscanf(fid_prot,'%g  ', [1 n_el_prot-1]); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    seg_index=0;
    %the total number of segments
    N_SEG = n_el_av_prot-1;
    
    %determine the length of a segment in pixel
    SEG_LENGTH = spline_p_av_prot(2) - spline_p_av_prot(1);
    
    %allocate memory
    av_score=zeros(N_SEG -START_SEG -END_SEG,1); 
    
    %figure
    %hold on
    for seg = START_SEG:N_SEG - END_SEG
        seg_index=seg_index+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% read segment boundary coordinates %%%%%%%%%%%%%%%
        if seg == START_SEG
            fileName_seg=char(filelist_seg(time));
            fid_bound_coord = fopen(fileName_seg,'r');       
            if fid_bound_coord == -1
                error('Could not open segment boundary coordinates');
            end  
            for i=1:START_SEG-1
                n_el_bound = fscanf(fid_bound_coord,'%g  ', [1 1]);   
                xb         = fscanf(fid_bound_coord,'%g  ', [1 n_el_bound]);
                yb         = fscanf(fid_bound_coord,'%g  ', [1 n_el_bound]);
            end          
        end 
        n_el_bound = fscanf(fid_bound_coord,'%g  ', [1 1]);   
        xb         = fscanf(fid_bound_coord,'%g  ', [1 n_el_bound]);
        yb         = fscanf(fid_bound_coord,'%g  ', [1 n_el_bound]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% create the same segment placed some way %%%%%%%%
        %%%%%%%%%%%%%%%%% set back along the normal %%%%%%%%%%%%%%%%%%%%%%
        xb_shift  = xb - SEG_SHIFT * x_av_normal(seg);
        yb_shift  = yb - SEG_SHIFT * y_av_normal(seg);  
        xb_shift2 = xb - SEG_SHIFT * x_av_normal(seg)- PROT_DEPTH * x_av_normal(seg);
        yb_shift2 = yb - SEG_SHIFT * y_av_normal(seg)- PROT_DEPTH * y_av_normal(seg);
        %plot(yb, xb,   '.b');
        %plot(yb_shift, xb_shift,   '.g');
        %plot(yb_shift2, xb_shift2, '.r');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% PROTRUSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        av_prot(seg_index)      = sqrt(x_av_prot(seg)^2 + y_av_prot(seg)^2);
        %protrusion into the normal direction
        scal                    = x_av_normal(seg) .* x_av_prot(seg) + y_av_normal(seg) .* y_av_prot(seg);
        x_av_prot_proj          = scal * x_av_normal(seg);
        y_av_prot_proj          = scal * y_av_normal(seg);
        av_prot_proj(seg_index) = sign(scal)*sqrt(x_av_prot_proj^2 + y_av_prot_proj^2);
        %get the angle between the normal direction and the retrograde flow
        av_prot_proj_angle(seg_index) = 180/pi * sign(scal) * real(acos(scal / av_prot(seg_index)));     
        %rescale!
        av_prot(seg_index)      = av_prot(seg_index)        * PIXEL / TIME_INTERVAL;
        av_prot_proj(seg_index) = av_prot_proj(seg_index)   * PIXEL / TIME_INTERVAL;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Object intensity values %%%%%%%%%%%%%%%%%%%%%%%%%  
        % extract the values
        img=zeros(imgSize);
        mask_image=roipoly(img,xb,yb);
        
        %extract from the image
        image_crop = mask_image .* image;
        
        clear image_crop_extr;
        [image_crop_extr(:,2), image_crop_extr(:,1), image_crop_extr(:,3)] = find(image_crop);
        
        %get number of measurments in segment
        av_cell_value_num(seg_index) = size(image_crop_extr,1);
        if av_cell_value_num(seg_index) > 0
            av_cell_value(seg_index) = sum(image_crop_extr(:,3));
            %rescale! This is score per second!
            av_cell_value(seg_index) = av_cell_value(seg_index) / av_cell_value_num(seg_index) / TIME_INTERVAL;
        else
            av_cell_value(seg_index)         = -99;
            av_cell_value_num(seg_index)     = -99; 
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         
         
        %%%%%%%%%%%%%%%%%%  Nice plot    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PROT_PLOT
            if seg == START_SEG
                h_prot_plot = figure;
                plot(x_p_edge, y_p_edge,'-r');
                hold on
                plot(x_p_edge_next,  y_p_edge_next, '-g');
                
%                 h_comp_plot = figure;
%                 plot(x_p_edge, y_p_edge,'-r');
%                 hold on  
%                 %get the edge described by the spline        
%                 max_par = (edge_sp_array_x(time_index).knots(end));
%                 p_n = 1:max_par;
%                 spline_pixel=[ fnval(edge_sp_array_y(time_index),p_n); fnval(edge_sp_array_x(time_index),p_n)]';
%                 plot(spline_pixel(:,2),spline_pixel(:,1),'.g','MarkerSize',0.5);
            end
            
            figure(h_prot_plot);
            hold on         
            quiver(x_av_pos_prot(seg), y_av_pos_prot(seg), x_av_prot(seg), y_av_prot(seg), 1,'r');    
            h_axes = gca; 
            set(h_axes,'xgrid','off','ygrid','off','box','off','Visible','off');
            axis equal
            axis([0 imgSize(2) 0 imgSize(1)])  
            axis ij
            
            if seg == N_SEG - END_SEG
                %hgsave(h_prot_plot,[dir_alpha 'prot_vec' num2str(time) '.fig']); 
                %print(h_prot_plot, [dir_alpha 'prot_vec' num2str(time) '.eps'],'-depsc2','-tiff'); 
                print(h_prot_plot, [dir_alpha 'prot_vec' num2str(time) '.tif'],'-dtiff'); 
                close(h_prot_plot);
            end           
        end
        
        
        
        
        if 0
            if time == FIRST_TIME+1 |  time == FIRST_TIME+7  %: 1: FIRST_TIME + TOTAL_TIME_STEPS -1
                if seg == START_SEG 
                    %h_score_img = figure; 
                    %imshow(speckle_img,[]);
                    %imwrite(speckle_img);
                    h_edge_nice = figure;
                    plot(x_p_edge, y_p_edge,'r');
                    h_axes = gca;  
                    set(h_axes,'xgrid','off','ygrid','off','box','off','Visible','off');
                    axis equal
                    axis([0 imgSize(2) 0 imgSize(1)])  
                    axis ij
                    
                    %h_score_nice = figure;
                    h_normals_nice  = figure;
                    h_retro_nice    = figure;
                    h_prot_nice     = figure;
                end

                %draw retrograde flow
                figure(h_retro_nice);
                if seg == START_SEG 
                    plot(x_p_edge, y_p_edge,'r'); 
                    hold on
                end
                %draw average protrusion
                figure(h_prot_nice);
                quiver(x_av_pos_prot(seg), y_av_pos_prot(seg), x_av_prot(seg), y_av_prot(seg), 40,'g');
                hold on
                h_axes = gca; 
                set(h_axes,'xgrid','off','ygrid','off','box','off','Visible','off');
                axis equal
                axis([0 imgSize(2) 0 imgSize(1)])  
                axis ij
                %draw average normal direction
                figure(h_normals_nice);
                quiver(x_av_pos_normal(seg), y_av_pos_normal(seg), x_av_normal(seg), y_av_normal(seg),25,'m');
                hold on
                h_axes = gca; 
                set(h_axes,'xgrid','off','ygrid','off','box','off','Visible','off');
                axis equal
                axis([0 imgSize(2) 0 imgSize(1)])  
                axis ij
                
                %save the nice plots
                if seg == N_SEG - END_SEG
                    print(h_normals_nice, [dir_alpha 'nice_normals.eps'],'-depsc2'); 
                    print(h_prot_nice,    [dir_alpha 'nice_prot.eps'],'-depsc2');
                    print(h_edge_nice,    [dir_alpha 'nice_edge.eps'],'-depsc2');  
                end
            end %plot data if
        end %if nice plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
          
    end %segment for
    
  
    if 0
        h_edge = figure;
        title('Edge position at last time step');
        plot(x_p_edge, y_p_edge,'r');       
    end

    if 0
        x_seg = START_SEG:1:N_SEG - END_SEG;
        figure
        plot(x_seg, alpha, 'r');
        hold on
        plot(x_seg, av_retrograde_flow_crop_proj, 'b:');
        plot(x_seg, av_prot_proj,'g-.');
        plot(x_seg, av_score,'b--');
        title(['Alpha at timestep  ' num2str(time)]);
        xlabel('Segment');
        legend('alpha','retrograde flow','protrusion','polymerization');
        
        figure
        plot(x_seg, score_sample_num,'b');
        hold on
        plot(x_seg, num_seg_velocities,'g');
        plot(x_seg, var_seg_velocities,'r');
        title(['Number of scores per segment at time  ' num2str(time)]); 
        text(imgSize(1),imgSize(2)-50,['Number of averaged time steps   '  num2str(time)]); 
        legend('number scores','number vel vectors','variance');
        xlabel('Segment');
    end
    
 
    %%%%%%%%%%%% Save the edge for the next time step
    n_pix       = n_pix_next;
    x_p_edge    = x_p_edge_next;
    y_p_edge    = y_p_edge_next;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Save it for time  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% The variables have the structure [segments,time] %%%%%%%%%%
    %%%%%%  These values are already RESCALED!!  %%%%%%%%%%%%%%%%%%%%%
    %protrusion
    protrusion_t(:,time_index)              = av_prot';
    %protrusion in normal direction
    protrusion_n_t(:,time_index)            = av_prot_proj';
    angle_prot_norm_t(:,time_index)         = av_prot_proj_angle';
    
    cell_value_t(:,time_index)              = av_cell_value';
    cell_value_num_t(:,time_index)          = av_cell_value_num';  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %time  for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Postprocessing  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  Remove first time step, because most probably  %%%%%%%%%%%
%%%%%%  there is no data                               %%%%%%%%%%%
protrusion_t(:,1)               = [];
protrusion_n_t(:,1)             = [];
angle_prot_norm_t(:,1)          = []; 

cell_value_t(:,1)            = [];  
cell_value_num_t(:,1)        = [];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Convert to images   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img_protrusion              = imresize(protrusion_n_t,            IMAGE_STRETCH, 'nearest');
img_cell_value              = imresize(cell_value_t,              IMAGE_STRETCH, 'nearest');
%%%%%%%%%%%   Filter              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_kernel = fspecial('gaussian',[IMG_GAUSS_W IMG_GAUSS_W], IMG_GAUSS_SIG);
img_protrusion_f            = imfilter(img_protrusion,      filter_kernel);
img_cell_value_f            = imfilter(img_cell_value,      filter_kernel);
%save these images
save([dir_alpha 'var_images.mat'], 'img_protrusion', 'img_protrusion_f', 'img_cell_value', 'img_cell_value_f');


h_protrusion_img = figure;
imshow(img_protrusion_f);
title('Protrusion');
xlabel('Time');
ylabel('Velocity [nm/s]');
colormap(jet);
if AUTO_COL_SCALE
    MAX_COLOR_VAL_PROT = max(img_protrusion_f(:));
end
caxis([-MAX_COLOR_VAL_PROT MAX_COLOR_VAL_PROT]);
hgsave(h_protrusion_img,[dir_alpha 'prot_activity_img.fig']);
print(h_protrusion_img, [dir_alpha 'prot_activity_img.eps'],'-depsc2','-tiff');
print(h_protrusion_img, [dir_alpha 'prot_activity_img.tif'],'-dtiff');

h_cell_value_img = figure;
imshow(img_cell_value_f);
title('Protein activity');
xlabel('Time');
ylabel('Velocity [nm/s]');
colormap(jet);
if AUTO_COL_SCALE
    MAX_COLOR_VAL_PROT = max(img_cell_value_f(:));
end
caxis([-MAX_COLOR_VAL_PROT MAX_COLOR_VAL_PROT]);
hgsave(h_cell_value_img,[dir_alpha 'prot_activity_img.fig']);
print(h_cell_value_img, [dir_alpha 'prot_activity_img.eps'],'-depsc2','-tiff');
print(h_cell_value_img, [dir_alpha 'prot_activity_img.tif'],'-dtiff');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot the variable space %%%%%%%%%%%%%%%%%%%%%
if VAR_SPACE_PLOT
    h_var_space = figure;
    hold on
%     c=colormap(colorcube(length(scores_t(:,1))));
    for t=1:TOTAL_TIME_STEPS - FIRST_TIME
         scatter3(scores_t(:,t), retrograde_flow_n_t(:,t), protrusion_n_t(:,t), 'filled');   
    end 
%    title('Color corresponds to segment');
    xlabel('Scores [s^{-1}]');
    ylabel('Retrograde Flow [nm/s]');
    zlabel('Protrusion [nm/s]');
    %save figure
    hgsave(h_var_space,[dir_alpha 'variable_space.fig']); 
    print(h_var_space, [dir_alpha 'variable_space.eps'],'-depsc2','-tiff'); 
    print(h_var_space, [dir_alpha 'variable_space.tif'],'-dtiff'); 
    hold off 
end %var_space_plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Save the variable to disk %%%%%%%%%%%%%%%%%%%
A=cat(3,protrusion_n_t);
save([dir_alpha 'A.mat'], 'A');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SURF_ACTIVITY_PLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%   Protrusion flow activity plot     %%%%%%%%%%%%%%
    max_val = max(protrusion_n_t(:));
    h_prot_activity = figure;
    temp = protrusion_n_t;
    surf(x_axis, y_axis, temp,...
        'EdgeLighting','phong','EdgeColor','none','FaceColor','interp');
    title('Protrusion');
    xlabel('Time [s]');
    ylabel('Edge position [\mu m]');
    zlabel('Protrusion [nm/s]]');  
    caxis([-max_val max_val]);
    colorbar;
    clear max_val
    clear temp
    %save figure
    hgsave(h_prot_activity,[dir_alpha 'prot_activity.fig']); 
    print(h_prot_activity, [dir_alpha 'prot_activity.eps'],'-depsc2','-tiff'); 
    print(h_prot_activity, [dir_alpha 'prot_activity.tif'],'-dtiff');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear x_axis;
    clear y_axis;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Bootstrap method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if BOOTSTRAP
    [bootstat_score_prot,bootsam]  = bootstrp(400,'corrcoef', scores_t,            protrusion_n_t);
    [bootstat_flow_prot, bootsam]  = bootstrp(400,'corrcoef', retrograde_flow_n_t, protrusion_n_t);
    [bootstat_flow_score, bootsam] = bootstrp(400,'corrcoef', retrograde_flow_n_t, scores_t);
    figure
    title('Correlation test with bootstrap method');
    subplot(3,1,1),hist(bootstat_score_prot(:,2));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w');
    title('Score - protrusion');
    hold on
    subplot(3,1,2),hist(bootstat_flow_prot(:,2));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w');
    title('Retrograde flow - protrusion');
    subplot(3,1,3),hist(bootstat_flow_score(:,2));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w');
    title('Retrograde flow - score');
end %bootstrap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% The variables have the structure [segments,time] %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Calculate TIME averaged values %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Good for cell with strong variations in segments %%%%%%%%%%
seg_samples = sum(protrusion_t ~= 0, 2);
av_protrusion_seg              = sum(protrusion_t,2)                ./ seg_samples;
av_protrusion_n_seg            = sum(protrusion_n_t,2)              ./ seg_samples;
av_angle_prot_norm_seg         = sum(angle_prot_norm_t,2)           ./ seg_samples;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Calculate SEGMENT averaged values %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Good for cell with strong variations in time %%%%%%%%%%%%%%
clear index;
index = protrusion_t ~= 0;
time_samples = sum(index,1);
av_protrusion_time              = sum(protrusion_t,1)               ./ time_samples;
av_protrusion_n_time            = sum(protrusion_n_t,1)             ./ time_samples;
av_angle_prot_norm_time         = sum(angle_prot_norm_t,1)          ./ time_samples;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Calculate total averages (scalar values) %%%%%%%%%%
av_total_prot_n             = mean(av_protrusion_n_seg);
av_total_prot               = mean(av_protrusion_seg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the axes used for ploting the data
x_seg_axis  = START_SEG:N_SEG - END_SEG;
x_time_axis = (FIRST_TIME: 1: TOTAL_TIME_STEPS - 1).* TIME_INTERVAL;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot TIME average values %%%%%%%%%%%%%%%%%%%%%%%%
h_time_av = figure;
subplot(3,1,1);
plot(x_seg_axis, av_protrusion_n_seg,'-og');
hold on
plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_prot_n,          '-g');
text(START_SEG+1, 1.3*av_total_prot_n,        num2str(av_total_prot_n));
title(['Segment statistics']);
ylabel('Velocity [nm/s]');
legend('Protrusion');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot the angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_angle = figure;
plot(x_seg_axis, av_angle_prot_norm_seg, '-db') 
xlabel('Segment #');
ylabel('Angle [degree]');
legend('Protrusion - normal direction');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear x_seg_axis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot SEGMENT average values %%%%%%%%%%%%%%%%%%%%%
h_seg_av = figure;
subplot(3,1,1);plot(x_time_axis, av_protrusion_n_time,'-og');
hold on
plot(x_time_axis, ones(size(x_time_axis)).*av_total_prot_n,                  '-g');
text(min(x_time_axis), 1.4*av_total_prot_n,             num2str(av_total_prot_n));
title(['Time statistics']);
ylabel('Velocity [nm/s]');
legend('Protrusion');

set(gcf,'Position', [10 10 560 700]);
hgsave(h_seg_av,[dir_alpha 'segment_av_variables.fig']); 
print(h_seg_av, [dir_alpha 'segment_av_variables.eps'],'-depsc2','-tiff'); 
print(h_seg_av, [dir_alpha 'segment_av_variables.tif'],'-dtiff');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot the angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_angle = figure;
plot(x_time_axis, av_angle_prot_norm_time, '-db') 
xlabel('Image #');
ylabel('Angle [degree]');
legend('Protrusion - normal direction');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear x_time_axis;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    Plot comparison of normal direction values %%%%%%%%%%
%%%%%%%%%    and non normal direction values            %%%%%%%%%%
figure
x_axis = START_SEG:N_SEG - END_SEG;
plot(x_axis, av_protrusion_n_seg,'-ob');
hold on
plot(x_axis, av_protrusion_seg,'-xb');
line([START_SEG length(av_protrusion_seg)], [av_total_prot_n av_total_prot_n],'Color','b');
line([START_SEG length(av_protrusion_seg)], [av_total_prot av_total_prot],'Color','b');
text(START_SEG+1, 1.2*av_total_prot_n,        num2str(av_total_prot_n));
text(START_SEG+1, 1.2*av_total_prot,          num2str(av_total_prot));
title(['Comparison between normal and non normal values']);
legend('n protrusion', 'protrusion', 'n retrograde flow', 'retrograde flow');
xlabel('Segment #');
ylabel('Velocity [pixel]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ALL_DATA_CALC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%% Reshape data into a single time segment vector!  %%%%%%%
    protrusion_n_t          = reshape(protrusion_n_t, prod(size(protrusion_n_t)),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%% Calculate statistic over All data!  %%%%%%%%%%%%%
    mean_protrusion         = mean(protrusion_n_t);
    var_protrusion          = var(protrusion_n_t);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%% Fillter segment time data  %%%%%%%%%%%%%%%%%%%%%%
    F_WINDOW=10;
    F_SIGMA=0.6;
    w=gausswin(F_WINDOW,F_SIGMA);
    protrusion_t_ff         =filtfilt(w/length(w),1,protrusion_n_t); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    figure
    plot(protrusion_n_t);
    hold on
    plot(protrusion_t_ff,'r');
    xlabel('Time*segments');
    ylabel('Protrusion');
    title(['Protrusion mean=  ' num2str(mean_protrusion) '   Var=  ' num2str(var_protrusion)]);
end %calculations over all data points

fclose(fid_pixel_edge);
fclose(fid_av_normal);
fclose(fid_bound_coord);
fclose(fid_av_prot);
