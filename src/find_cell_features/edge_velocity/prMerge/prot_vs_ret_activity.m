function prot_vs_ret_stat_activity(edge_parameters, merg_parameters, post_parameters, var1, var2, var1_name, var2_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Compare protrusion vs retraction                              %
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
ACT_TIME_SHIFT      = post_parameters.act_time_shift;  
SCO_TIME_SHIFT      = post_parameters.delay;  


% Create result directory for correlation data
mkdir([RESULT_DIR filesep 'prot_vs_ret' filesep]);
RESULT_DIR = [RESULT_DIR filesep 'prot_vs_ret' filesep];

result_file_name = [var1_name '_' var2_name];

if merg_parameters.units == 0
    xLabelText = 'frame #';
    yLabelText = 'Velocity (pixel/frame)';
elseif merg_parameters.units == 1
    xLabelText = 'Time (s)';    
    yLabelText = 'Velocity (nm/s)';
else  merg_parameters.units == 2
    xLabelText = 'Time (min)';    
    yLabelText = 'Velocity (um/min)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  First get the NaN locations in the data            %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%var1_nan = isnan(reshape(var1,1,[]));
%var2_nan = isnan(reshape(var2,1,[]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filter data                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SCO_TIME_SHIFT > 0
    var2 = circshift(var2,[0 SCO_TIME_SHIFT]);
    var2(:,1:SCO_TIME_SHIFT)=[];
    var1(:,1:SCO_TIME_SHIFT)=[];
elseif SCO_TIME_SHIFT < 0
    var2 = circshift(var2,[0 SCO_TIME_SHIFT]);
    var2(:,end+SCO_TIME_SHIFT:end)=[];
    var1(:,end+SCO_TIME_SHIFT:end)=[];
end


var1_array = reshape(var1,1,[]);
var2_array = reshape(var2,1,[]);

% remove zeros
index_z = find(~var2_array);
var1_array(index_z) = [];
var2_array(index_z) = [];

% delete Nan's from the sequence
var1_nan = isnan(var1_array);
var2_nan = isnan(var2_array);
not_nan_index = find((var1_nan & var2_nan));
var1_array(not_nan_index) = [];
var2_array(not_nan_index) = [];

% get protrusion only
var1_array_p_i = find(var1_array>0);
var1_array_p = var1_array(var1_array_p_i);
var2_array_p = var2_array(var1_array_p_i);

% calculate the average time in protrusion
tpsum = sum(var1>0,2)* edge_parameters.time_interval;
AveTimeInProtrusion = mean(tpsum(tpsum~=0))/60;

% get retraction
var1_array_n_i = find(var1_array<0);
var1_array_n = var1(var1_array_n_i);
var2_array_n = var2_array(var1_array_n_i);

% calculate the average time in retraction
tnsum = sum(var1<0,2)* edge_parameters.time_interval;
AveTimeInRetraction = mean(tnsum(tnsum~=0))/60;

%save([RESULT_DIR 'Ave-time-in-Prot-Ret.mat'],'AveTimeInProtrusion','AveTimeInRetraction');%,'-ascii');
fid = fopen([RESULT_DIR filesep 'Ave_Times_in_Prot_Ret'] ,'w+');
fprintf(fid, 'AveTimeInProtrusion (min.) = %f\n',       AveTimeInProtrusion);
fprintf(fid, 'AveTimeInRetraction (min.) = %f',       AveTimeInRetraction);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Statistical analysis        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean
mean_var2_p = mean(var2_array_p);
mean_var2_n = mean(var2_array_n);
% std
std_var2_p = std(var2_array_p);
std_var2_n = std(var2_array_n);

% histogram
[n_var2_p, xout_var2_p] = hist(var2_array_p,40);
[n_var2_n, xout_var2_n] = hist(var2_array_n,40);

% t-test
[h_t, significance, ci] = ttest2(var2_array_p,var2_array_n,0.01,'both','unequal');
% ranksum test
[p,h_ranksum]           = ranksum(var2_array_p,var2_array_n,'alpha',0.01);

figure
plot(xout_var2_p, n_var2_p, 'r');
hold on
plot(xout_var2_n, n_var2_n, 'b');
line([mean_var2_p mean_var2_p],[0 max(n_var2_p)],'Color','r');
line([mean_var2_n mean_var2_n],[0 max(n_var2_n)],'Color','b');
title(['Test for equal means; t test:  ', num2str(h_t+0), ' Rank sum test:  ', num2str(h_ranksum+0)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% plot(var1_array_p,var2_array_p,'r+');
% hold on
% plot(var1_array_n,var2_array_n,'bo');
% 
% xlabel(yLabelText);
% ylabel('Activity [-]');

% the following plot identifies various segments at different
% location and time (x,t) that represent protrusion(pos-prot) or retraction
% (neg-protrusion) vs. the corresponding value of activity at those
% segments. Mohsen 2/6/08.

h_prot_vs_ret = figure;

plot(var1_array_p,var2_array_p,'r+');
hold on
plot(var1_array_n,var2_array_n,'bo');

xlabel(yLabelText);
% ylabel('Activity [-]');
ylabel(var2_name);
    
title('statistics of activity values correspoinding with protrusion or retraction');
legend(['pos-' var1_name 'vs. corresp ' var2_name 'values'], ['neg-' var1_name 'vs. corresp ' var2_name 'values'], 'Best');
    
hgsave(h_prot_vs_ret,[RESULT_DIR 'prot_ret.fig']);
print(h_prot_vs_ret, [RESULT_DIR 'prot_ret.eps'],'-depsc2','-tiff');
print(h_prot_vs_ret, [RESULT_DIR 'prot_ret.tif'],'-dtiff');    
    


