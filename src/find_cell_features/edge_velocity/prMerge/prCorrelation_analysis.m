function prCorrelation_analysis(edge_parameters, merg_parameters, post_parameters, var1, var2, var1_name, var2_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Estimate the auto correlations and the                        %
%           cross correlations                                            %
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
corr_lag            = post_parameters.corr_lag;
RESULT_DIR          = post_parameters.result_dir;

% Create result directory for correlation data
mkdir([RESULT_DIR filesep 'correlation' filesep]);
RESULT_DIR = [RESULT_DIR filesep 'correlation' filesep];


if merg_parameters.units == 0
    xLabelText = 'frame #';
    yLabelText = '(pixel frame) units';
elseif merg_parameters.units == 1
    xLabelText = 'Time (s)';    
    yLabelText = '(nm s) units';
else  merg_parameters.units == 2
    xLabelText = 'Time (min)';    
    yLabelText = '(um min) units';
end


PLOT_FIGURES = 1;

%if PEAK_VARIABILITY_ANALYSIS = -1 then it will do bootstrap(lines142-183)
%the default is "0" 
% Note: for this option make sure you uncomment the line for "return" 
      % statement at before the end of  if PEAK_VARIABILITY_ANALYSIS
      % operations. Also for default "0" option make sure the statements
      % "if 0" and the corresponding "end" are uncommented.
      
  PEAK_VARIABILITY_ANALYSIS = 0;


ALIGN_CORR         = -1 ;
%    method = 0 : cross correlation based alignment
%    method = 1 : maximum value based alignment
%    method = 2 : minimum value based alignment


% get the number of segments
% get the number of time steps
n_segments = size(var1, 1);
n_time     = size(var1, 2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filter data                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD > 1
    [var1_filtered, var1_noise] = prFilterTimeSeries(var1,FILTER_METHOD, windowSize, spline_s);
    [var2_filtered, var2_noise] = prFilterTimeSeries(var2,FILTER_METHOD, windowSize, spline_s);
    if trend == 1
            var1 = var1_filtered;
            var2 = var2_filtered;
    else
            var1 = var1_noise;
            var2 = var2_noise;       
    end
end

% Calculate the significance bounds
% 99% probability -> 2.58
% 95% probability -> 1.96
% 90% probability -> 1.645
upper_sig_bound = (n_time)^(-0.5) * 1.96;
lower_sig_bound = -upper_sig_bound;
upper_sig_bound_spatial = (n_segments)^(-0.5) * 1.96;
lower_sig_bound_spatial = -upper_sig_bound_spatial;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the correlation matrix and take the zero lag element the
% 'coeff', to normalize the sequence so the auto-correlation at zero lag
% are identically 1.0
% lag            = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lag  = corr_lag;
x_lag = -lag:lag;
time_interval = edge_parameters.time_interval;
x_lag = x_lag .* time_interval;
x_axis = time_interval.*(0:n_time-1);

dT = edge_parameters.time_interval;

for i_seg=1:n_segments
    corr_var1_var2(i_seg,:) = xcov(var1(i_seg,:), var2(i_seg,:), lag,'coeff');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot the average variotion of the activities, protrusion in time %%%%
%%%%% This might be useful for persistent motion %%%%%%%%%%%%%%%%%%%%%%%%%%
av_var1 = mean(var1);
av_var2 = mean(var2);
figure;
plot(x_axis, av_var1);
title(['Temporal variation of avarage '  var1_name ]);
xlabel('Time #');
ylabel('Average Value');

figure;
plot(x_axis, av_var2);
title(['Temporal variation of avarage '  var2_name ]);
xlabel('Time #');
ylabel('Average Value');

h_av_var1 = figure;
waterfall(x_axis, 1:size(corr_var1_var2,1), var1);
xlabel('Time');
ylabel('Segment');
zlabel(['average value of ' var1_name]);
hgsave(h_av_var1,   [RESULT_DIR 'av_' var1_name '_segments']);
% print(h_av_var1,    [RESULT_DIR 'av_' var1_name '_segments'],'-depsc2','-tiff');
% print(h_av_var1,    [RESULT_DIR 'av_' var1_name '_segments'],'-dtiff');

h_av_var2 = figure;
waterfall(x_axis, 1:size(corr_var1_var2,1), var2);
xlabel('Time');
ylabel('Segment');
zlabel(['average value of ' var2_name]);
hgsave(h_av_var2,   [RESULT_DIR 'av_' var2_name '_segments']);
% print(h_av_var2,    [RESULT_DIR 'av_' var2_name '_segments'],'-depsc2','-tiff');
% print(h_av_var2,    [RESULT_DIR 'av_' var2_name '_segments'],'-dtiff');


% for i_seg=1:n_segments
%     av_seg_var1(i_seg)=mean(var1(i_seg,:));
%     std_seg_var1(i_seg)=std(var1(i_seg,:));
%     av_seg_var2(i_seg)=mean(var2(i_seg,:));
%     std_seg_var2(i_seg)=std(var2(i_seg,:));
% end

av_seg_var1=mean(var1,2);
std_seg_var1=std(var1,0,2);
av_seg_var2=mean(var2,2);
std_seg_var2=std(var2,0,2);

clear var1_segment_statistics;
clear var2_segment_statistics;

var1_segment_statistics=[av_seg_var1, std_seg_var1];
save([RESULT_DIR 'av_seg_' var1_name], 'var1_segment_statistics');
var2_segment_statistics=[av_seg_var2, std_seg_var2];
save([RESULT_DIR 'av_seg_' var2_name], 'var2_segment_statistics');
figure;
errorbar(1:n_segments, av_seg_var1, std_seg_var1,'o');
xlabel('Segment');
ylabel(['Average ' var1_name ' in each segment'])
figure;
errorbar(1:n_segments, av_seg_var2, std_seg_var2,'o');
xlabel('Segment');
ylabel(['Average ' var2_name ' in each segment'])
%  save([RESULT_DIR 'av_seg_' var1_name ], '-ascii', 'av_seg_var1', 'std_seg_var1');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the following if statement chooses the alignment method defined in lines 41-43
%method 1 and 2 give the 2-D correlation map plus cross-correlation
%coefficients between edge and activity for each sampling window, i.e. as
%depicted in figures 2-K and 2-O of the multiplex paper. Method 0 gives the
%default 1-D cross-correlation of edge with activity.

if 0
    % Test output showing the time series and correlation of one window
    test_window = 30;
    %if PLOT_FIGURES
        figure
        plot(x_axis,var1(test_window,:),'g');
        hold on
        plot(x_axis,var2(test_window,:),'b');
        title(['Test plot '  var1_name '  ' var2_name ' in window ', num2str(test_window)]);
        xlabel('Time #');
        ylabel('Velocity');

        figure
        plot(x_lag, corr_var1_var2(test_window,:));
        title(['Test plot, correlation '  var1_name '  ' var2_name ' in window ', num2str(test_window)]);

        h_av_cov_seg = figure;
        plot(corr_var1_var2(:,lag+1));
        legend(var1_name, var2_name);
        title('Correlation between leading edge variables');
        xlabel('Segment');
        ylabel('Correlation');
   
        % plot correlation from all segments
        figure,waterfall(x_lag, 1:size(corr_var1_var2,1), corr_var1_var2);
        xlabel('Time shift');
        ylabel('Segment');
        zlabel('Correlation');
        
        %y_contour = 1:n_segments;
        %x_contour = x_lag;
        
        % the following 5 lines create the 2-D correlation map it doesn't
        % provide any new info. beyond the 1D cross correlations. Its just
        % for visualization and takes time, collaborators asked to comment it out.
        
%         [x_contour,y_contour] = meshgrid(x_lag,1:n_segments);
%         figure,contourf(x_contour,y_contour,corr_var1_var2,20);
%         xlabel('Time shift');
%         ylabel('Segment');    
%         colorbar;
        
     %end
end


if PEAK_VARIABILITY_ANALYSIS
    % interpolate splines into  each cross correlations
    for i_seg=1:n_segments
        corr_spline_interpolation{i_seg} = spapi(4,x_lag,corr_var1_var2(i_seg,:));
    end
    
    % plot the splines
    figure;
    hold on
    for i_seg=1:n_segments
        fnplt(corr_spline_interpolation{i_seg});
    end
    xlabel('Time shift');
    title('4th order spline interpolations');
    
    % Obtain peak from spline interpolations and error in peak estimation 
    % from the bootstrap
    for i_seg=1:n_segments
        % compute first derivative
        fprime = fnder(corr_spline_interpolation{i_seg});
        %fun = @(x)fprime;
        % find zeros
        time_lag(i_seg) = fzero(@(x) fun(x,fprime), 0);
        
    end
    figure
    plot(time_lag);
    disp('Mean time lag: ');
    mean(time_lag)
    disp('Std time lag: ');
    std(time_lag)
    
    
    % Try alternative method:
    % Fit smoothing spline to all data points
    % Get a estimation of the cross correlation curve
    y_data = reshape(corr_var1_var2', 1, prod(size(corr_var1_var2)));
    x_data = repmat(x_lag, 1, size(corr_var1_var2,1));
%     [max_x_residuals] = fit_GTPases(x_data, y_data);
    [max_x_residuals] = regression_bootstrap(x_data, y_data);

    return;

end


av_corr_var1_var2 = mean(corr_var1_var2,1);
% compute standard deviation of correlation 
se_corr_var1_var2 = std(corr_var1_var2,0,1);%/sqrt(size(corr_var1_var2,1));

h_av_cross_corr = figure;
plot(x_lag, av_corr_var1_var2)
hold on
plot(x_lag, av_corr_var1_var2+se_corr_var1_var2,'--');
plot(x_lag, av_corr_var1_var2-se_corr_var1_var2,'--');

line([x_lag(1); x_lag(end)],[upper_sig_bound; upper_sig_bound],'LineStyle','--');
line([x_lag(1); x_lag(end)],[lower_sig_bound; lower_sig_bound],'LineStyle','--');


%%%%% Note %%%%%%%%%
%the following 4 lines  are for putting the average cross-correlation 
%of each individual vindow on the 1-D cross-correlation diagram as in figure 2L
%of multiplexing paper. only the line that has the command "plot" is
%needed. The other 3 lines make a nicer plot with a gray cloud background.
% In fact line with bins = 25 (line 205) is irrelevant.
%%%% End Note %%%%%%%

%bins = 25;
%plot(repmat(x_lag,size(corr_var1_var2,1),1) , corr_var1_var2, 'bx');
%CMAP = flipud(gray(256));
%scattercloud(repmat(x_lag,size(corr_var1_var2,1),1) ,corr_var1_var2, 60, 1, 'k.', CMAP);

legend([var1_name ' - ' var2_name]);
title_str_1 = ['Averaged Cross correlation of ' var1_name ' and ' var2_name];
title([title_str_1]);
xlim([x_lag(1) x_lag(end)]);
xlabel('Lag (sec)');
ylabel('Correlation');
hgsave(h_av_cross_corr,   [RESULT_DIR 'cross_corr_' var1_name '_' var2_name]);
print(h_av_cross_corr,    [RESULT_DIR 'cross_corr_' var1_name '_' var2_name],'-depsc2','-tiff');
print(h_av_cross_corr,    [RESULT_DIR 'cross_corr_' var1_name '_' var2_name],'-dtiff');

% comment the following "return" out if you want the 1-D cross correlation and
% other plots that were usually generated in default choise of
% PEAK_VARIABILITY = 0.

%return

clear c;
c = [x_lag; av_corr_var1_var2];
save([RESULT_DIR 'cross_corr_' var1_name '_' var2_name '.mat'], 'c');

disp(['Number of time steps =', num2str(n_time)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                   %%%
%%%%%   Calculate the AUTO COVARIANCE to detect possible frequencies    %%%
%%%%%                                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_seg=1:n_segments
    auto_corr_var1(i_seg,:) = xcov(var1(i_seg,:), lag,'coeff');
    auto_corr_var2(i_seg,:) = xcov(var2(i_seg,:), lag,'coeff');    
end

% Average over segments
av_auto_corr_var1 = mean(auto_corr_var1,1);
av_auto_corr_var2 = mean(auto_corr_var2,1);

h_av_auto_corr = figure;
plot(x_lag, av_auto_corr_var1,'g')
hold on
plot(x_lag, av_auto_corr_var2,'r')
line([x_lag(1); x_lag(end)],[upper_sig_bound; upper_sig_bound],'LineStyle','--');
line([x_lag(1); x_lag(end)],[lower_sig_bound; lower_sig_bound],'LineStyle','--');
xlim([x_lag(1) x_lag(end)]);
legend(var1_name, var2_name);
title(['Averaged auto-correlation of ' var1_name ' and ' var2_name]);
xlabel('Time Lag ');
ylabel('Auto-correlation coefficient');
hgsave(h_av_auto_corr,   [RESULT_DIR 'auto_corr_' var1_name '_' var2_name]);
print(h_av_auto_corr,    [RESULT_DIR 'auto_corr_' var1_name '_' var2_name],'-depsc2','-tiff');
print(h_av_auto_corr,    [RESULT_DIR 'auto_corr_' var1_name '_' var2_name],'-dtiff');

%xlswrite([RESULT_DIR 'auto_corr_' var1_name '_' var2_name], [x_lag; av_auto_corr_var1; av_auto_corr_var2], 'cross_correlation');
c = [x_lag; av_auto_corr_var1; av_auto_corr_var2];
save([RESULT_DIR 'auto-correlation_' var1_name '_' var2_name '.mat'],'c');

clear c;
c = [x_lag; av_auto_corr_var1];
save([RESULT_DIR 'auto-correlation_' var1_name '.mat'],'c');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calculate spatial auto-correlation in var1 and var2  %%%%%%%%%%%%%%
%%%%%%% together with spatial cross-correlation btw var1 & var2 %%%%%%%%%%%

% s_lag = -10:10; % Doesn't work, has to be same size as the half of no. segments
%slag = round(n_segments/2);
slag = round(n_segments);
s_lag = -slag:slag;
% s_lag =round(-n_segments/2):round(n_segments/2);
ss_lag=s_lag.*post_parameters.segment_length_av./1000;


for i_time=1:n_time
    spatial_auto_corr_var1(i_time,:) = xcov(var1(:,i_time), slag,'coeff');
    spatial_auto_corr_var2(i_time,:) = xcov(var2(:,i_time), slag,'coeff');
    spatial_cross_corr_var1_var2(i_time,:) = xcov(var1(:,i_time), var2(:,i_time), slag, 'coeff');
% % spatial_cross_corr_var1_var2 is the spatial cross correlation between the variable 1 and
% variable 2.

    
    sp1(i_time) = csaps(s_lag, spatial_auto_corr_var1(i_time,:), 1);
    fn1(i_time) = fncmb(sp1(i_time),'-',0.2);
    tmp1 = fnzeros(fn1(i_time));
    if size(tmp1,2) > 0
        l_c1(i_time) = tmp1(1,end/2+1);
    else
        l_c1(i_time) = 0;
    end
    
    sp2(i_time) = csaps(s_lag, spatial_auto_corr_var2(i_time,:), 1);
    fn2(i_time) = fncmb(sp2(i_time),'-',0.2);
    tmp2 = fnzeros(fn2(i_time));
    if size(tmp2,2) > 0
        l_c2(i_time) = tmp2(1,end/2+1);
    else
        l_c2(i_time) = 0;
    end    
end


% The spatial correlation length development
figure
plot(l_c1*post_parameters.segment_length_av./1000,'g');
hold on
plot(l_c2*post_parameters.segment_length_av./1000,'r');
title('Characteristic length development (microns)');
xlabel('Frame #');
ylabel('0.2 correlation');
legend(var1_name, var2_name);


 [x_contour,y_contour] = meshgrid(ss_lag,(1:n_time)*time_interval);
%         figure,surf(x_contour,y_contour,spatial_auto_corr_var1);
%         title(['Time evoloution of spatial Auto correlation of ' var1_name]);
%         xlabel('Spatial shift (Micron)');
%         ylabel('Time (s)');    
%         colorbar;
%         figure,surf(x_contour,y_contour,spatial_auto_corr_var2);
%         title(['Time evoloution of spatial Auto correlation of ' var2_name]);
%         xlabel('spatial shift (Micron)');
%         ylabel('Time (s)');    
%         colorbar;
 h_spatial_cross_corr_var1_var2 = figure;
%         pcolor(x_contour,y_contour,spatial_cross_corr_var1_var2);
        surf(x_contour,y_contour,spatial_cross_corr_var1_var2);
        title(['Time evoloution of spatial Cross-correlation of ' var1_name ' - ' var2_name]);
        xlabel('Spatial shift (Micron)');
        ylabel('Time (s)');    
        colorbar;
 hgsave(h_spatial_cross_corr_var1_var2,   [RESULT_DIR 'spatial_temporal_cross_corr_' var1_name '_' var2_name '.fig']);
 print(h_spatial_cross_corr_var1_var2,    [RESULT_DIR 'spatial_temporal_cross_corr_' var1_name '_' var2_name '.eps'],'-depsc2','-tiff');
 print(h_spatial_cross_corr_var1_var2,    [RESULT_DIR 'spatial_temporal_cross_corr_' var1_name '_' var2_name '.tif'],'-dtiff');
%         figure,contourf(x_contour,y_contour,spatial_cross_corr_var1_var2);
%         xlabel('spatial shift (Micron)');
%         ylabel('Time (s)');    
%         colorbar;
%

% Average over time
av_spatial_auto_corr_var1 = mean(spatial_auto_corr_var1,1);
av_spatial_auto_corr_var2 = mean(spatial_auto_corr_var2,1);
av_spatial_cross_corr_var1_var2 = mean(spatial_cross_corr_var1_var2,1);
% compute standard deviation of correlation
se_spatial_cross_corr_var1_var2 = std(spatial_cross_corr_var1_var2,0,1);%/sqrt(size(corr_var1_var2,1));


h_av_spatial_cross_corr_var1_var2 = figure;
plot(ss_lag, av_spatial_cross_corr_var1_var2);
hold on;
plot(ss_lag, av_spatial_cross_corr_var1_var2+se_spatial_cross_corr_var1_var2,'--');
plot(ss_lag, av_spatial_cross_corr_var1_var2-se_spatial_cross_corr_var1_var2,'--');
line([ss_lag(1); ss_lag(end)],[upper_sig_bound_spatial; upper_sig_bound_spatial],'LineStyle','--');
line([ss_lag(1); ss_lag(end)],[lower_sig_bound_spatial; lower_sig_bound_spatial],'LineStyle','--');

title(['Average spatial Cross correlation of ' var1_name ' - ' var2_name]);
legend([ var1_name ' - ' var2_name ]);
xlabel('spatial shift (Micron)');
ylabel('Cross-correlation coefficient');
hgsave(h_av_spatial_cross_corr_var1_var2,   [RESULT_DIR 'spatial_cross_corr_' var1_name '_' var2_name '.fig']);
print(h_av_spatial_cross_corr_var1_var2,    [RESULT_DIR 'spatial_cross_corr_' var1_name '_' var2_name '.eps'],'-depsc2','-tiff');
print(h_av_spatial_cross_corr_var1_var2,    [RESULT_DIR 'spatial_cross_corr_' var1_name '_' var2_name '.tif'],'-dtiff');

h_av_spatial_auto_corr = figure;
plot(ss_lag, av_spatial_auto_corr_var1, 'g');
hold on
plot(ss_lag, av_spatial_auto_corr_var2, 'r');
title(['Average spatial Auto correlation of ' var1_name ' and ' var2_name]);
legend(var1_name, var2_name);
xlabel('Spatial shift (Micron)');
ylabel('Auto-correlation coefficient');
hgsave(h_av_spatial_auto_corr,   [RESULT_DIR 'spatial_auto_corr_' var1_name '_' var2_name '.fig']);
print(h_av_spatial_auto_corr,    [RESULT_DIR 'spatial_auto_corr_' var1_name '_' var2_name '.eps'],'-depsc2','-tiff');
print(h_av_spatial_auto_corr,    [RESULT_DIR 'spatial_auto_corr_' var1_name '_' var2_name '.tif'], '-dtiff');

clear c;
c = [ss_lag; av_spatial_cross_corr_var1_var2];
save([RESULT_DIR 'spatial_cross_corr_' var1_name '_' var2_name '.mat'], 'c');

clear c;
c = [ss_lag; av_spatial_auto_corr_var1; av_spatial_auto_corr_var2];
save([RESULT_DIR 'spatial_auto-correlation_' var1_name '_' var2_name '.mat'], 'c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function f = paraOpt(p, x, y)
% p = [a,b,m]

f = y - (p(1)*(x - p(3)).^2 + p(2));

function f  = fun(x, spline)
f = fnval(x,spline);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


