function visualize_activity(edge_parameters, merg_parameters, post_parameters, var1, var2, var1_name, var2_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Visualize the data                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MERG_DIR            = post_parameters.merg_dir;  
PROJECT_DIR         = post_parameters.project_dir;  
POST_PROC_OPTIONS   = post_parameters.post_proc_options;   
DO_FILTER           = post_parameters.do_filter;
FILTER_METHOD       = post_parameters.filter_method;
trend               = post_parameters.trend;
windowSize          = post_parameters.windowSize;
spline_s            = post_parameters.spline_s;
RESULT_DIR          = post_parameters.result_dir;
% function specific
IMAGE_STRETCH       = post_parameters.image_stretch;
ISOTRPIC            = post_parameters.isotropic;  
IMG_GAUSS_W         = post_parameters.img_gauss_w;  
IMG_GAUSS_SIG       = post_parameters.img_gauss_sig;  
TIME_INTERVAL       = edge_parameters.time_interval;

UNITS = merg_parameters.units;
% % determine what units to use

if UNITS == 0
    xLabelText = 'frame #';
    yLabelText = 'Velocity (pixel/frame)';
elseif UNITS == 1
    xLabelText = 'Time (s)';    
    yLabelText = 'Velocity (nm/s)';
else  UNITS == 2
    xLabelText = 'Time (min)';    
    yLabelText = 'Velocity (um/min)';
end
% x_time_axis = (merg_parameters.first_time: 1: merg_parameters.total_time_steps).* TIME_INTERVAL;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filter data                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD > 1
    [var1_filterd, var1_noise] = prFilterTimeSeries(var1, FILTER_METHOD, windowSize, spline_s);
    [var2_filterd, var2_noise] = prFilterTimeSeries(var2, FILTER_METHOD, windowSize, spline_s);

    if trend == 1
        var1 = var1_filterd;
        var2 = var2_filterd;
    else
        var1 = var1_noise;
        var2 = var2_noise;           
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

img_var1  = imresize(var1, IMAGE_STRETCH, 'nearest');
img_var2  = imresize(var2, IMAGE_STRETCH, 'nearest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%            Filter              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FILTER = 1;
if FILTER
    isotropic_filter_kernel     = fspecial('gaussian',[IMG_GAUSS_W IMG_GAUSS_W], IMG_GAUSS_SIG);
    anisotropic_filter_kernel   = anisotropicGaussD(IMG_GAUSS_SIG/10, IMG_GAUSS_SIG, 0, 3);

    if ISOTRPIC
        filter_kernel = isotropic_filter_kernel;
    else
        filter_kernel = anisotropic_filter_kernel;
    end

    [img_var1, dum]    = Gauss2DBorder(img_var1, IMG_GAUSS_SIG);
    [img_var2, dum]    = Gauss2DBorder(img_var2, IMG_GAUSS_SIG);
end

  
%%%%% modified by Mohsen Sabouri-Ghomi on 7/19/07 for putting time stamp %%%%%%%%
%%%%         and Segment # stamp on the activity maps             %%%%%%%%%

%%%% set the no. of ticks to be put on the activity map %%%%%%%%%
No_Ticks = 10;
iptsetpref('ImshowAxesVisible','on');

h_var1 = figure;
imshow(img_var1,[]);
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
colormap(jet);
colorbar;
title([var1_name '  '  yLabelText]);
xlabel('Time (s)');
ylabel('Segment #');

%hgsave(h_var1,[RESULT_DIR 'activity_' var1_name '.fig']);
%print(h_var1, [RESULT_DIR 'activity_' var1_name '.eps'],'-depsc2','-tiff');
%print(h_var1, [RESULT_DIR 'activity_' var1_name '.tif'],'-dtiff');

% Get the zero level
zero_level = contourc(img_var1, [0 0]);


h_var2 = figure;
imshow(img_var2,[]);
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
colormap(hot);
colorbar;
hold on
plot(zero_level(1,:),zero_level(2,:),'w.','MarkerSize',3);
title([var2_name ' (Intensity/Pixel)']);
xlabel('Time (s)');
ylabel('Segment #');

%hgsave(h_var2,[RESULT_DIR 'activity_' var2_name '.fig']);
%print(h_var2, [RESULT_DIR 'activity_' var2_name '.eps'],'-depsc2','-tiff');
%print(h_var2, [RESULT_DIR 'activity_' var2_name '.tif'],'-dtiff');

figure
surface(img_var1);
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
xlabel('Time (s)');
ylabel('Segment #');
zlabel(yLabelText);