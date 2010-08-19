function prStatistical_analysis(edge_parameters, merg_parameters, post_parameters, var1, var2, var1_name, var2_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyze the data                                              %
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

% Create result directory for correlation data
mkdir([RESULT_DIR filesep 'statistical_analysis' filesep]);
RESULT_DIR = [RESULT_DIR filesep 'statistical_analysis' filesep];


if merg_parameters.units == 0
    xLabelText = 'frame #';
    yLabelText = 'Velocity (pixel/frame)';
    TimeUnit = 1/edge_parameters.time_interval;
elseif merg_parameters.units == 1
    xLabelText = 'Time (s)';    
    yLabelText = 'Velocity (nm/s)';
    TimeUnit = 1;
else  merg_parameters.units == 2
    xLabelText = 'Time (min)';    
    yLabelText = 'Velocity (um/min)';
    TimeUnit = 1/60;
end

%x_time_axis = (merg_parameters.first_time: 1: merg_parameters.total_time_steps).* edge_parameters.time_interval;
x_time_axis = getTimeAxis(merg_parameters,edge_parameters.time_interval,TimeUnit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%            Detrend   data                          %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD > 1
    [var1_filtered, var1_noise] = prFilterTimeSeries(var1,FILTER_METHOD, windowSize, spline_s);
    [var2_filtered, var2_noise] = prFilterTimeSeries(var2,FILTER_METHOD, windowSize, spline_s);
else
    var1_filtered = var1;
    var2_filtered = var2;
    var1_noise = var1;
    var2_noise = var2;
end
    
av_time_var1_filtered = mean(var1_filtered,1);
av_time_var1_noise    = mean(var1_noise,1);
av_time_var1          = mean(var1,1);

av_time_var2_filtered = mean(var2_filtered,1);
av_time_var2_noise    = mean(var2_noise,1);
av_time_var2          = mean(var2,1);


h_var1 = figure;
plot(x_time_axis, av_time_var1_filtered,'g');
hold on
plot(x_time_axis, av_time_var1,'g:');
xlim([x_time_axis(1) x_time_axis(end)]);
title(['Filtered segment averaged ' var1_name]);
xlabel(xLabelText);
hgsave(h_var1,   [RESULT_DIR ['seg_av_' var1_name] '.fig']);
print(h_var1,    [RESULT_DIR ['seg_av_' var1_name] '.eps'],'-depsc2','-tiff');
print(h_var1,    [RESULT_DIR ['seg_av_' var1_name] '.tif'],'-dtiff');

% calculate rms of filtered protrusion_n
% this is modified by Mohsen to generate Root Mean Square of normal
% protrusion velocity in each segment.

if (strcmp(var1_name,'protrusion_n')==1)
    START_SEG       = merg_parameters.start_seg;
    START_SEG = START_SEG+1;
    END_SEG         = merg_parameters.end_seg;
    SEG_NR          = merg_parameters.seg_nr;
    x_seg_axis  = START_SEG:SEG_NR - END_SEG;
    for i=1:size(var1_filtered,1)
        rms_var1_filtered_seg(i)    = norm(var1_filtered(i,:))/sqrt(size(var1_filtered,2));
    end
    h_rms_var1_filtered = figure;
    plot(x_seg_axis, rms_var1_filtered_seg,'g');
    legend(['RMS of ' var1_name], 'Location','Best');
    title('Root Mean Square variation of Protrusion_n in each segment');
    xlabel('Segment #');
    ylabel(yLabelText);
    hgsave(h_rms_var1_filtered,   [RESULT_DIR ['rms_' var1_name] '.fig']);
    print(h_rms_var1_filtered,    [RESULT_DIR ['rms_' var1_name] '.eps'],'-depsc2','-tiff');
    print(h_rms_var1_filtered,    [RESULT_DIR ['rms_' var1_name] '.tif'],'-dtiff');
end
%end calculate rms of filtered protrusion_n

h_var2 = figure;
plot(x_time_axis, av_time_var2_filtered,'g');
hold on
plot(x_time_axis, av_time_var2,'g:');
xlim([x_time_axis(1) x_time_axis(end)]);
title(['Filtered segment averaged ' var2_name]);
xlabel(xLabelText);
hgsave(h_var2,   [RESULT_DIR ['seg_av_' var2_name] '.fig']);
print(h_var2,    [RESULT_DIR ['seg_av_' var2_name] '.eps'],'-depsc2','-tiff');
print(h_var2,    [RESULT_DIR ['seg_av_' var2_name] '.tif'],'-dtiff');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot normalized mean curves
h_var1_var2 = figure;
plot(x_time_axis, prestd(av_time_var1_filtered),'g');
hold on
plot(x_time_axis, prestd(av_time_var2_filtered),'r');
title('Averaged curves, normalized to mu=0, sigma=1');
legend(var1_name, var2_name);
xlabel(xLabelText);
hgsave(h_var1_var2,   [RESULT_DIR ['normalized_' var1_name ' - ' var2_name] '.fig']);
print(h_var1_var2,    [RESULT_DIR ['normalized_' var1_name ' - ' var2_name] '.eps'],'-depsc2','-tiff');
print(h_var1_var2,    [RESULT_DIR ['normalized_' var1_name ' - ' var2_name] '.tif'],'-dtiff');

% % plot curves of window 12
% figure;
% plot(x_time_axis, var1_filtered(12,:),'g');
% title('Curves window 12');
% legend(var1_name);
% xlabel(xLabelText);
% 
% figure;
% plot(x_time_axis, var2_filtered(12,:),'r');
% % title('Curves window 12 normalized to mu=0, sigma=1');
% title('Curves window 12');
% legend(var2_name);
% xlabel(xLabelText);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var1_array = reshape(var1,prod(size(var1)),1); 
var2_array = reshape(var2,prod(size(var2)),1); 

var1_std = std(var1_array);
var2_std = std(var2_array);

var1_mean = mean(var1_array);
var2_mean = mean(var2_array);

[var1_hist,xout] = hist(var1_array,50);
h_var1_hist = figure;    
plot(xout, var1_hist);
hold on
line([var1_mean, var1_mean], [0 max(var1_hist)],'LineStyle',':');
line([var1_mean - var1_std, var1_mean - var1_std], [0 max(var1_hist)],'LineStyle',':');
line([var1_mean + var1_std, var1_mean + var1_std], [0 max(var1_hist)],'LineStyle',':');
title([var1_name '  histogram', '  Width: ', num2str(2*var1_std)]);

[var2_hist,xout] = hist(var2_array,50);
h_var2_hist = figure;    
plot(xout, var2_hist);
hold on
line([var2_mean, var2_mean], [0 max(var2_hist)],'LineStyle',':');
line([var2_mean - var2_std, var2_mean - var2_std], [0 max(var2_hist)],'LineStyle',':');
line([var2_mean + var2_std, var2_mean + var2_std], [0 max(var2_hist)],'LineStyle',':');
title([var2_name '  histogram', '  Width: ', num2str(2*var2_std)]);

% This is to define data distribution, outliers etc. with boxplot
%figure;
%boxplot of protrusion
%boxplot(var1_array);
%title(['Boxplot of ' var1_name]);
figure;
%boxplot of activity
boxplot(var2_array);
title(['Boxplot of ' var2_name]);

[H, pValue] = lillietest(var1_array,0.05) 
[H, pValue] = lillietest(var2_array,0.05) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Normalized histogramms
h_area =  trapz(xout,var1_hist);
var1_hist = var1_hist./ h_area; 

h_area =  trapz(xout,var2_hist);
h_norm_hist = figure;
var2_hist = var2_hist./ h_area; 
plot(xout, var1_hist);
hold on
plot(xout, var2_hist,'r');
legend(var1_name, var2_name);
title('Normalized  histogram');
hgsave(h_norm_hist,   [RESULT_DIR ['normalized_histo_' var1_name ' - ' var2_name] '.fig']);
print(h_norm_hist,    [RESULT_DIR ['normalized_histo_' var1_name ' - ' var2_name] '.eps'],'-depsc2','-tiff');
print(h_norm_hist,    [RESULT_DIR ['normalized_histo_' var1_name ' - ' var2_name] '.tif'],'-dtiff');



function x_time_axis = getTimeAxis(merg_parameters,time_interval,TimeUnit)

if isfield(merg_parameters,'protTimePoints')
   x_time_axis = merg_parameters.protTimePoints.*time_interval.*TimeUnit;
elseif isfield(merg_parameters,'flowTimePoints')
   x_time_axis = merg_parameters.flowTimePoints.*time_interval.*TimeUnit;
elseif isfield(merg_parameters,'scoreTimePoints')
   x_time_axis = merg_parameters.scoreTimePoints.*time_interval.*TimeUnit;
elseif isfield(merg_parameters,'activityTimePoints')
   x_time_axis = merg_parameters.activityTimePoints.*time_interval.*TimeUnit;
elseif isfield(merg_parameters,'activity_TimePoints')
   x_time_axis = merg_parameters.activity_TimePoints.*time_interval.*TimeUnit;
else
   x_time_axis = (merg_parameters.first_time+1: 1: merg_parameters.total_time_steps).*time_interval;
end
