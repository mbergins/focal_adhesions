function prMergeProtrusion(merg_parameters, prot_parameters)
% PRMERGE MERGES PROTRUSION
%
%
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
% Matthias Machacek 08/04/06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% units:
% 0: pixel/frame
% 1: pixel/frame
% 2: um/min
UNITS           = merg_parameters.units;
PIXEL           = prot_parameters.pixel;
TIME_INTERVAL   = prot_parameters.time_interval;
PROJ_FLOW       = 1;
PROT_SAMPLING   = prot_parameters.prot_sampling;
PARENTH_R       = prot_parameters.parenth_r;
PARENTH_L       = prot_parameters.parenth_l;
PROJECT_DIR     = merg_parameters.project_dir;
PROT_DIR        = merg_parameters.prot_dir;
PLOT_VECTOR_MODE= merg_parameters.plot_vector_mode;
SEG_NR          = merg_parameters.seg_nr;
TOTAL_TIME_STEPS= merg_parameters.total_time_steps;
START_SEG       = merg_parameters.start_seg;
START_SEG = START_SEG+1;
END_SEG         = merg_parameters.end_seg;
PLOT_VECTORS    = merg_parameters.plot_vectors;



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
iEntry = 1;
for( i = 1:length(filelist_temp))
   if(~filelist_temp(i).isdir)
      filelist_mask(iEntry) = {[PROJECT_DIR filesep PROT_DIR filesep 'cell_mask' filesep filelist_temp(i).name]};
      iEntry = iEntry + 1;
   end
end
clear filelist_temp


% If we generate overlay images we have to determine which images to use
if PLOT_VECTORS
    % get image file names
    [IMG_DIR, name, ext, versn] = fileparts(prot_parameters.file);
    
    %[IMG_DIR, name, ext, versn] = fileparts('S:\scripps\analysis\machacek\FocalAdhesions\80505\set_2\pax2_cut1\crop_default001.tif');
    
    filelist_images_b = dir([IMG_DIR filesep '*.tif']);
    if length(filelist_images_b) == 0
        filelist_images_b = dir([IMG_DIR filesep '*.bmp']);
    end

    iEntry = 1;
    for i=1:length(filelist_images_b)
        if(~filelist_images_b(i).isdir)
            filelist_images(iEntry,:) = filelist_images_b(i).name;
            iEntry = iEntry + 1;
        end
    end
end

% file name of the the pixel edge
file_pixel_edge = [PROJECT_DIR  filesep PROT_DIR filesep 'pixel_edge.mat'];

% file name of the normal vectors
file_normal=[PROJECT_DIR filesep PROT_DIR filesep 'normal_matrix.mat'];

% get the filelist filenames containing the protrusion vectors    
file_protrusion=[PROJECT_DIR filesep PROT_DIR filesep 'protrusion.mat'];
%%%%%%%%%%%%%% End get file name list  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% load data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(file_pixel_edge);

% load data with edge splines
load([PROJECT_DIR filesep  PROT_DIR filesep 'edge_spline.mat']);

% this loads cell normal_matrix 
load(file_normal);

% this loads cell 'protrusion'
load(file_protrusion);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%the first time step refers to number of the files in the folder% 
%and NOT to the number in the file name!!!!!!!!!!!!   %%%%%%%%%%%
time_index = 0;
strg = sprintf('%%.%dd',3);
FIRST_TIME = 1;
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
    %%%%%  average the unit normals in segments %%%%%%%%%%%%%%%%%%%%%%
    % average the unit normals in segments
    p_n=1:3:n_pix; % the number 3 comes from the imEdgeTracker.m
    [x_normal_out_av, y_normal_out_av, x_av_pos_normal, y_av_pos_normal, m_pos]=...
        prGetAvEdge(edge_sp_array_x(time), edge_sp_array_y(time), edge_sp_array_x(time).knots(end),...
        p_n, normal_matrix{time}(:,3), normal_matrix{time}(:,4), 'nr_sect', SEG_NR);

    % re-normalize the averaged normals
    l_n_av = sqrt(x_normal_out_av.^2 + y_normal_out_av.^2);
    x_av_normal = x_normal_out_av ./ l_n_av;
    y_av_normal = y_normal_out_av ./ l_n_av;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  average the protrusion in segments %%%%%%%%%%%%%%%%%%%%%%%%
    
    % Average the protrusion vectors in the segments
    i_nn=1:PROT_SAMPLING:n_pix;
    l_i_nn=length(i_nn);
    i_nn(l_i_nn-PARENTH_R:l_i_nn)=[];
    i_nn(1:PARENTH_L)=[];
    clear l_i_nn;

    [x_av_prot, y_av_prot, x_av_pos_prot, y_av_pos_prot, spline_p_av_prot]=...
        prGetAvEdge(edge_sp_array_x(time), edge_sp_array_y(time), n_pix,...
        i_nn, protrusion{time}(:,3), protrusion{time}(:,4), 'nr_sect', SEG_NR);

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
        [xb, yb] = prCreateSegments(edge_sp_array_x(time), edge_sp_array_y(time),...
                     x_av_normal, y_av_normal, imgSize(2), imgSize(1), SEG_NR, 2, seg);        
        %%%%%%%%%%%% End Generate segments  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% PROTRUSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%  Nice plot    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PLOT_VECTORS
            if PLOT_VECTOR_MODE == 1
                % all in one plot
                if seg == START_SEG && seg_index == 1
                    h_nice = figure;
                    %overlay_image = imread(overlay_image_path);
                    overlay_image = imread([IMG_DIR filesep filelist_images(time_index,:)]);
                    imshow(overlay_image,[0 2300]);
                    hold on
                    plot(x_p_edge, y_p_edge,'r');
                    h_axes = gca;
                    set(h_axes,'xgrid','off','ygrid','off','box','off','Visible','off');
                    axis equal
                    axis([0 imgSize(2) 0 imgSize(1)])
                    axis ij
                end
                figure(h_nice);
                plot(xb,yb,'y','LineWidth',0.5);
                plot(xb_shift,  yb_shift,'y','LineWidth',0.5);
                if seg == START_SEG
                    text(xb(1),yb(1),'1','Color','w');
                elseif seg == SEG_NR - END_SEG
                    text(xb(end),yb(end), num2str(SEG_NR - END_SEG),'Color','w'); 
                end

                quiver(x_av_pos_prot(seg), y_av_pos_prot(seg), x_av_prot(seg), y_av_prot(seg), 10,'g');
             
                if seg == SEG_NR - END_SEG
                    hh = getframe(gcf);
                    imwrite(hh.cdata,[merg_dir 'overlay' filesep 'overlay_' filelist_images(time_index,:)],'tif');
                    close(h_nice);
                end
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


clear img;
clear mask_image;
clear mask_img;
clear velocity_field;
clear velocity_field_crop;
clear x_p_edge;
clear y_p_edge;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Postprocessing  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Re-scale segment length
segment_length = segment_length .* PIXEL;
% Determine average segment length
segment_length_av = mean(segment_length);
figure
plot(segment_length)
title('Segment length (nm)');
save([merg_dir 'segment_length_av'],'segment_length_av');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Save the variable to disk %%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([merg_dir 'protrusion.mat'], 'protrusion_seg');
save([merg_dir 'protrusion_n.mat'], 'protrusion_normal');
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
av_protrusion_seg               = sum(protrusion_seg,2)   ./ sum(protrusion_seg ~=0, 2);
av_protrusion_n_seg             = sum(protrusion_normal,2)./ sum(protrusion_normal ~=0, 2);
av_protrusion_angle_seg         = sum(protrusion_angle,2) ./ sum(protrusion_angle ~=0, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Calculate SEGMENT averaged values %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Good for cell with strong variations in time %%%%%%%%%%%%%%
av_protrusion_time              = sum(protrusion_seg,1)      ./ sum(protrusion_seg ~=0, 1);
av_protrusion_n_time            = sum(protrusion_normal,1)   ./ sum(protrusion_normal ~=0, 1);
av_protrusion_angle_time        = sum(protrusion_angle,1)    ./ sum(protrusion_angle ~=0, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Calculate total averages (scalar values) %%%%%%%%%%%%%%%%%%%
av_total_prot_n             = mean(av_protrusion_n_seg);
av_total_prot               = mean(av_protrusion_seg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the x-axis used for ploting the data
x_seg_axis  = START_SEG:SEG_NR - END_SEG;
x_time_axis = (FIRST_TIME: 1: TOTAL_TIME_STEPS).* TIME_INTERVAL;


if (TOTAL_TIME_STEPS - FIRST_TIME) > 5
    %%%%%%%%%%%%%%   Calculate a robust estimage of protrusion %%%%%%%
    [robust_av_protrusion, stats] = robustfit(x_time_axis, av_protrusion_n_time);
    robust_sigma_protrusion = stats.robust_s;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %no robust regression can be performed if there are too few
    %parameters!
    robust_av_protrusion(1) =               av_total_prot_n;
    robust_av_protrusion(2) =               0;
    robust_sigma_protrusion =               0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Save the average variables to disk %%%%%%%%%%

save([av_data_dir 'av_prot_variables.mat'], 'av_protrusion_n_seg','av_protrusion_n_time',...
     'robust_av_protrusion', 'robust_sigma_protrusion');

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
elseif UNITS == 1
    xLabelText = 'Time (s)';    
    yLabelText = 'Velocity (nm/s)';
else
    xLabelText = 'Time (min)';    
    yLabelText = 'Velocity (um/min)';
end

h_red_time_av = figure;
plot(x_seg_axis, av_protrusion_n_seg,'g');
hold on
plot(x_seg_axis, ones(size(x_seg_axis)).*av_total_prot_n,'-g');
text(min(x_seg_axis), 1.4*av_total_prot_n,             num2str(av_total_prot_n));
xlim([x_seg_axis(1) x_seg_axis(end)]);
legend('Protrusion');
title(['Time statistics']);
xlabel('Segment #');
ylabel(yLabelText);

hgsave(h_red_time_av,[figures_dir 'time_av_prot_retro.fig']);
print(h_red_time_av, [figures_dir 'time_av_prot_retro.eps'],'-depsc2','-tiff');
print(h_red_time_av, [figures_dir 'time_av_prot_retro.tif'],'-dtiff');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot the angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_angle = figure;

plot(x_seg_axis, av_protrusion_angle_seg, '-db')
xlabel('Segment #');
ylabel('Angle [degree]');
legend( 'Protrusion - normal direction');

hgsave(h_angle,[figures_dir 'angle.fig']);
print(h_angle, [figures_dir 'angle.eps'],'-depsc2','-tiff');
print(h_angle, [figures_dir 'angle.tif'],'-dtiff');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plot SEGMENT average values %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_red_seg_av = figure;
plot(x_time_axis, av_protrusion_n_time,'g');
hold on
legend('Protrusion');
plot(x_time_axis, ones(size(x_time_axis)).*av_total_prot_n,'-g');
text(min(x_time_axis), 1.4*av_total_prot_n, num2str(av_total_prot_n));
title(['Segment statistics']);
xlabel(xLabelText);
ylabel(yLabelText);


hgsave(h_red_seg_av,[figures_dir 'segment_av_prot_retro.fig']);
print(h_red_seg_av, [figures_dir 'segment_av_prot_retro.eps'],'-depsc2','-tiff');
print(h_red_seg_av, [figures_dir 'segment_av_prot_retro.tif'],'-dtiff');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plot the angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_angle = figure;

plot(x_time_axis, av_protrusion_angle_time, '-db')
xlabel('Image #');
ylabel('Angle [degree]');
legend('Protrusion - normal direction');

hgsave(h_angle,[figures_dir 'angle_time_av.fig']);
print(h_angle, [figures_dir 'angle_time_av.eps'],'-depsc2','-tiff');
print(h_angle, [figures_dir 'angle_time_av.tif'],'-dtiff');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    Plot comparison of normal direction values %%%%%%%%%%
%%%%%%%%%    and non normal direction values            %%%%%%%%%%
x_axis = START_SEG:SEG_NR - END_SEG;

figure
plot(x_axis, av_protrusion_n_seg,'-ob');
hold on
plot(x_axis, av_protrusion_seg,'-xb');
line([START_SEG length(av_protrusion_seg)], [av_total_prot_n av_total_prot_n],'Color','b');
line([START_SEG length(av_protrusion_seg)], [av_total_prot av_total_prot],'Color','b');
text(START_SEG+1, 1.2*av_total_prot_n,        num2str(av_total_prot_n));
text(START_SEG+1, 1.2*av_total_prot,          num2str(av_total_prot));
title(['Comparison between normal and non normal values']);
legend('n protrusion', 'protrusion');
xlabel('Segment #');
ylabel('Velocity [pixel]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Convert to images   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (TOTAL_TIME_STEPS - FIRST_TIME) > 5
    img_protrusion  = imresize(protrusion_normal,IMAGE_STRETCH, 'nearest');

    [img_protrusion_f, dum] = Gauss2DBorder(img_protrusion, IMG_GAUSS_SIG);
    save([merg_dir 'img_prot.mat'], 'img_protrusion',  'img_protrusion_f');

    h_protrusion_img = figure;
    imshow(img_protrusion_f);
    title('Protrusion');
    xlabel(xLabelText);
    ylabel('Cell edge');
    colormap(jet);
    MAX_COLOR_VAL_PROT = abs(robust_av_protrusion(1)+0.5*robust_av_protrusion(2)*x_time_axis(end)) + 3*robust_sigma_protrusion;
    caxis([-MAX_COLOR_VAL_PROT MAX_COLOR_VAL_PROT]);
    colorbar;
    hgsave(h_protrusion_img,[figures_dir 'prot_activity_img.fig']);
    print(h_protrusion_img, [figures_dir 'prot_activity_img.eps'],'-depsc2','-tiff');
    print(h_protrusion_img, [figures_dir 'prot_activity_img.tif'],'-dtiff');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose('all');








