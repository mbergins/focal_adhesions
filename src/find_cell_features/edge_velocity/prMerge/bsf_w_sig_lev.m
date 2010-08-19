function bsf(string, n, var1_name, var2_name, nimages)
% fits a smoothing spline into correlation coefficient obtained from varous
% cells and provides inferential statistics such as Confidence intervals, and
% location of the extrema (maximum or minimum) correlation together with  
% standard error of the extremum location.
% 
%
%
% SYNOPSIS      bsf('common_file_name', n, 'var1_name', 'var2_name', nimages)
%
% INPUT                :  'commonFilename' 
%                      : example : 'corr_protrusion_n_activity1_'
%                           
%                      :  n ---> this is the number of files with common 
%                         name representing the cross-correlation of two 
%                         activities from "n" cells.
%                      : 'var1_name' and 'var2_name' are the names of the 
%                         two activities that you have correlated for the n
%                         individual cells.
%                      example: 'Protrusion', 'Cdc42' etc.
%                      : 'nimages' is the no. of images you entered in
%                         PrMergePanel GUI when you merged your data.
%               Note I :  no. of images for all cells should be the same for
%                         this to give you correct upper and lower 
%                         significance levels. 
%                        
%              Note II :  if the no. of segments (windows) along the edge 
%                         for each individual cell you calculated the spatial
%                         correlation for were not equal to each other, you
%                         will get an error. There is a work around for
%                         this which requires trimming your spatial
%                         correlations to match the correlation of the cell
%                         with minimum length. 
%              
%                      
% OUTPUT               : 
%
% Mohsen Sabouri Ghomi 8/12/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear global f_hat;
clear global smoothing_factor;
clear global CrossCorr_var1_var2;
clear global x_lag;
clear global n;
clear global x_data;
clear global y_data;
clear global nimages;
clear global upper_sig_bound;
clear global lower_sig_bound;

global n;
global f_hat;
global x_data;
global y_data;
global smoothing_factor;
global CrossCorr_var1_var2;
global x_lag;
global upper_sig_bound; 
global lower_sig_bound;
global nimages;





load ([string num2str(1)]);
 lng=size(c,2);
 clear c;
 c=zeros(2,lng);


 for i=1:n
     load ([string num2str(i)]);
     CrossCorr_var1_var2(i,:) = c(2,:);
     %     clear c
 end
 x_lag = c(1,:);

 y_data = reshape(CrossCorr_var1_var2', 1, prod(size(CrossCorr_var1_var2)));
 x_data = repmat(x_lag, 1, size(CrossCorr_var1_var2,1));
 [extrema_x_residuals] = regression_bootstrap(x_data, y_data,var1_name,var2_name);

 
end %bsf function;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% regression_bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%

function [extrema_x_residuals] = regression_bootstrap(x_data, y_data,var1_name,var2_name)
%
global n;
global f_hat;
global x_data;
global y_data;
global smoothing_factor;
global CrossCorr_var1_var2;
global x_lag;
global upper_sig_bound; 
global lower_sig_bound;
global nimages;


% put data into right shape
if size(x_data,1) < size(x_data,2)
   x_data =  x_data';
end
if size(y_data,1) < size(y_data,2)
   y_data =  y_data';
end


% bootstrapping pairs (x_i, y_i)
%BOOTSTRAP_MODE = 1;
% bootstrapping residuals r
%BOOTSTRAP_MODE = 2;

BOOTSTRAP_MODE = 2; 

% Number of bootstrap samples 

nboot = 2000;

% Confidence level: alpha
alpha = 0.05;

% spline smoothing factor: p 
% p=1/(1+h^3/0.6). h=ave(data spacing). here h=xx_step.
% GTPases setting : smoothing_factor = 0.0001;
smoothing_factor = 0.85;
% smoothing_factor = 0.5;

% Calculate the significance bounds
% 99% probability -> 2.58
% 95% probability -> 1.96
% 90% probability -> 1.645
% av_n_time = 184; 
% this the average time series for Src-K treated cells
%av_n_time = 103.3; 
% this is the average time series for control cells
% this the average number of time steps (image frames) of each cell n.
%upper_sig_bound = (av_n_time)^(-0.5) * 1.96;
%lower_sig_bound = -upper_sig_bound;
upper_sig_bound = (nimages)^(-0.5) * 1.96;
lower_sig_bound = -upper_sig_bound;

xx_min = min(x_data);
xx_max = max(x_data);
xx_step = (xx_max - sign(xx_min)*abs(xx_min))/100;
if xx_step > 1
    xx_step = 1;
end
xx = xx_min:xx_step:xx_max;

h_indi_spline = figure;
waterfall(x_lag, 1:n,CrossCorr_var1_var2);
xlabel('Time shifts (sec)');
ylabel('cell #');
zlabel('Correlation coefficient');
title(['individual ' var1_name ' & ' var2_name ' correlations together with the fit spline (shown at cell # 0)']);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a estimation of the curve f: f_hat
f_hat = csaps(x_data, y_data, smoothing_factor);
% Get a estimation of the curve f: f_hat, with a Loess algorithm
span = 5;
% y_hat_loess = smooth(x_data, y_data, 0.1 ,'rloess');
% y_hat_sgolay = smooth(x_data, y_data,'sgolay');

% compute estimated predictor values at regular x positions
y_hat = ppval(f_hat,xx);
waterfall(xx, 0, ppval(f_hat, xx));
hgsave(h_indi_spline, 'correaltions-fitspline');
hold off;

% plot the smoothing spline with data points
% figure
% plot(xx, ppval(f_hat, xx));
% hold on
% plot(x_data, y_data, 'x');
% %plot(x_data, y_hat_loess, ':');
% %plot(x_data, y_hat_sgolay, ':');
% title('Data with spline fit')'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now calculate the residuals using the estimation
residuals = y_data - ppval(f_hat, x_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%         Analyse residuals         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Get rms
rms = norm(residuals)/sqrt(length(residuals));
disp(['SD of residuals =', num2str(rms)]);

figure
plot(x_data, residuals, 'x');
title('Residuals');

h_hist = figure;
hist(residuals,60);
title('Histogram of residuals');
% Estimate the density function of the residuals
[res_density, xi] = ksdensity(residuals);
hgsave(h_hist, 'residuals_histogram');

h_pdf = figure;
plot(xi, res_density);
title('Estimate of probability density funct. (pdf) from residuals');
hgsave(h_pdf, 'pdf_residuals');
print(h_pdf,  'pdf_residuals','-dpng');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      Bootstrap to compute the statistics     %%%%%%%%%%%%%%%%
if BOOTSTRAP_MODE ==1
    % Bootstrap pairs (x,y) to calculate a f_star
    bootstat = bootstrp(nboot, @bootfun1, x_data, y_data);
    
elseif BOOTSTRAP_MODE == 2 
    % Use the residuals for bootstraping and to calculate a f_star by 
    % computing a new y_star: y_star = f_star(x_data) + bootstrapped(residual)
    bootstat = bootstrp(nboot, @bootfun2, residuals);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Compute all bootstrap regression estimates   %%%%%%%%%
% I commented the following two lines as it was plottin every booted sample
% and therefore crowding the fig. 5. Instead I added these tow lines after
% the following end. Mohsen Sabouri

h_boot_samples = figure;  
hold on 
for i=1:nboot
    y_star(i,:) = ppval(bootstat(i), xx);
    % Extract position and value of peak for each f_star
    [max_y_star(i), max_index] = max(y_star(i,:));
    max_x_star(i) = xx(max_index);
    %the following line is for geting the minimum of y_star
%     [min_y_star(i), min_x_star] = fnmin(bootstat(i));
    [min_y_star(i), min_index] = min(y_star(i,:));
    min_x_star(i) = xx(min_index);
    plot(xx,y_star(i,:));
end
% figure %by M
% hold on %by M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%       Compare with estimate of f: f_hat

plot(xx, y_hat, 'r');
title('Bootstrap samples of correlation coefficients together with spline fit');
hgsave(h_boot_samples, ['bootstrap_samples_' var1_name '_' var2_name]);
print(h_boot_samples, ['bootstrap_samples_' var1_name '_' var2_name], '-dpng');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute residual at all xx locations
for i=1:nboot
    y_residuals(i,:) = y_hat - y_star(i,:);
end
    
% Compute Standard deviation for each xx location
std_star = std(y_star,0,1);
% Compute Standard error for each xx location
se_star = std_star./sqrt(nboot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate alpha/2 and 1-alpha/2 confidence intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The percentile interval
y_alpha_l = quantile(y_residuals, alpha/2,   1);
y_alpha_u = quantile(y_residuals, 1-alpha/2, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BCA method for confidence intervals
% p 184
 if 0
% bias correaction z_0
for i=1:length(xx)
    z_0 = norminv(sum(y_star(:,i) < y_hat(i))/nboot)

    alpha1 = 0;
    alpha2 = 0;
end

% acceleration a
for i=1:length(x_data)
    x_jn = x_data;
    y_jn = y_data;    
    x_jn(i) = [];
    y_jn(i) = [];    
    
    % Compute f_hat_i
    f_hat_i = csaps(x_jn, y_jn, smoothing_factor);    
    
end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%       Analyse maxima    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get peak of f_hat
% the following 2 lines were used originally
[max_y_hat, max_index] = max(y_hat);
max_x_hat = xx(max_index);

% [max_y_hat, max_x_hat] =  fnmin(fncmb(f_hat,-1));
% max_y_hat = -max_y_hat;
[min_y_hat, min_x_hat] = fnmin(f_hat);

maxyhat=max_y_hat
minyhat=min_y_hat
absmax_y_hat = abs(max_y_hat)
absmin_y_hat = abs(min_y_hat)

 if (abs(max_y_hat) >= abs(min_y_hat))
%  if (abs(max_y_hat) < abs(min_y_hat))   

    % compute residuals of the x and y positions of the peak
    max_x_residuals = max_x_hat - max_x_star;
    max_y_residuals = max_y_hat - max_y_star;
    extrema_x_residuals = max_x_residuals;

    % compute std
    std_max_x_star = std(max_x_star,0);
    std_max_y_star = std(max_y_star,0);

    % Compute Standard error for each xx location
    se_max_x_star = std_max_x_star./sqrt(nboot);
    se_max_y_star = std_max_y_star./sqrt(nboot);

    % alpha/2  and 1-alpha/2  confidence interval
    max_x_alpha_l = quantile(max_x_residuals, alpha/2);
    max_x_alpha_u = quantile(max_x_residuals, 1-alpha/2);

    max_y_alpha_l = quantile(max_y_residuals, alpha/2);
    max_y_alpha_u = quantile(max_y_residuals, 1-alpha/2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Plot results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot standard error pluss spline fit to the data
h_boot = figure;
hold on;
plot(x_data, y_data, 'gd', 'LineWidth',1.3);% 'MarkerFaceColor',[.49 1 .63]);
plot(xx, y_hat, 'r'); 
plot(xx, y_hat + se_star, 'c--');
plot(xx, y_hat - se_star, 'c--');

line([xx(1); xx(end)],[upper_sig_bound; upper_sig_bound],'LineStyle','--');
line([xx(1); xx(end)],[lower_sig_bound; lower_sig_bound],'LineStyle','--');

% title('spline fit to data from cells together with standard errors from bootstraping')
% hold off;
% plot confidence interval pluss spline fit to the data
% figure
% % plot(xx, y_hat, 'b');
% hold on;
plot(xx, y_hat + y_alpha_u, 'r--');
plot(xx, y_hat + y_alpha_l, 'r--');
plot(max_x_hat, max_y_hat,'black');



% plot standard error
myErrorbar(max_x_hat, max_y_hat,[se_max_x_star se_max_y_star], [se_max_x_star  se_max_y_star]);

% plot confidence interval
myErrorbar(max_x_hat, max_y_hat,[max_x_alpha_l max_y_alpha_u], [max_x_alpha_l  max_y_alpha_u]);
% title('spline fit to data from cells together with 95% confidence intervals from bootstraping')
xlabel('Time shifts (sec)');
ylabel('Correlation coefficient');
title(['Bootstraped correlation coefficient between ' var1_name ' & ' var2_name] );
hgsave(h_boot,   ['bootstrap_' var1_name '_' var2_name]);
print(h_boot, ['bootstrap_' var1_name '_' var2_name], '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the individual data but colored
% plot(x_data(1:41), y_data(1:41), 'black');
% plot(x_data(1:41), y_data(42:82), 'blue');
% plot(x_data(1:41), y_data(83:123), 'green');
% plot(x_data(1:41), y_data(124:164), 'yellow');
% plot(x_data(1:41), y_data(165:205), 'magenta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' Peak location')
disp(['X  = ' , num2str(max_x_hat), ' ', num2str(max_x_alpha_l) , '  +',num2str(max_x_alpha_u)]);
disp(['Y  = ' , num2str(max_y_hat), ' ', num2str(max_y_alpha_l) , '  +',num2str(max_y_alpha_u)]);
disp(['Confidence interval= ', num2str(100*(1-alpha)), '%']);
disp(['Bootstrap samples= ', num2str(nboot)]);
disp(['Spline smoothing factor= ', num2str(smoothing_factor)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save x-position
save('max_x_star','max_x_star');

else

 %%%%%%%%%%%%%%%%%%%%      Analyse minima    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute residuals of the x and y positions of the minimum neg. corr.
    min_x_residuals = min_x_hat - min_x_star;
    min_y_residuals = min_y_hat - min_y_star;
    extrema_x_residuals = min_x_residuals;

    % compute std
    std_min_x_star = std(min_x_star,0);
    std_min_y_star = std(min_y_star,0);

    % Compute Standard error for each xx location
    se_min_x_star = std_min_x_star./sqrt(nboot);
    se_min_y_star = std_min_y_star./sqrt(nboot);

    % alpha/2  and 1-alpha/2  confidence interval
    min_x_alpha_l = quantile(min_x_residuals, alpha/2);
    min_x_alpha_u = quantile(min_x_residuals, 1-alpha/2);

    min_y_alpha_l = quantile(min_y_residuals, alpha/2);
    min_y_alpha_u = quantile(min_y_residuals, 1-alpha/2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Plot results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot standard error pluss spline fit to the data
h_boot = figure;
hold on;
plot(x_data, y_data, 'gd', 'LineWidth',1.3);% 'MarkerFaceColor',[.49 1 .63]);
plot(xx, y_hat, 'r'); 
plot(xx, y_hat + se_star, 'c--');
plot(xx, y_hat - se_star, 'c--');

line([xx(1); xx(end)],[upper_sig_bound; upper_sig_bound],'LineStyle','--');
line([xx(1); xx(end)],[lower_sig_bound; lower_sig_bound],'LineStyle','--');

% title('spline fit to data from cells together with standard errors from bootstraping')
% hold off;
% plot confidence interval pluss spline fit to the data
% figure
% plot(xx, y_hat, 'b');  % plots the spline fit
% hold on;
 plot(xx, y_hat + y_alpha_u, 'r--'); % upper confidence interval
 plot(xx, y_hat + y_alpha_l, 'r--'); % lower confidence interval
 plot(min_x_hat, min_y_hat,'black'); % shows the min point as one data point



% plot standard error
myErrorbar(min_x_hat, min_y_hat,[se_min_x_star se_min_y_star], [se_min_x_star  se_min_y_star]); % this is burried under the confidence interval

% plot confidence interval
% title('spline fit to data from cells together with 95% confidence intervals from bootstraping') original title
%myErrorbar(min_x_hat, min_y_hat,[min_x_alpha_l min_y_alpha_u],[min_x_alpha_l  min_y_alpha_u]); % this appears on top of standard error plot
%myErrorbar(min_x_hat, min_y_hat,[min_x_alpha_l min_y_alpha_l],[min_x_alpha_u  min_y_alpha_u]);% this plot the carculated lower and upper 95% conf. intervals (unsymmetric). 
myErrorbar(min_x_hat, min_y_hat,[min_x_alpha_l min_y_alpha_u]);%,[min_y_alpha_l  min_y_alpha_u]);


xlabel('Time shifts (sec)');
ylabel('Correlation coefficient');
title(['Bootstraped correlation coefficient between ' var1_name ' & ' var2_name] );
hgsave(h_boot,   ['bootstrap_' var1_name '_' var2_name]);
print(h_boot, ['bootstrap_' var1_name '_' var2_name], '-dpng');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the individual data but colored
% plot(x_data(1:41), y_data(1:41), 'black');
% plot(x_data(1:41), y_data(42:82), 'blue');
% plot(x_data(1:41), y_data(83:123), 'green');
% plot(x_data(1:41), y_data(124:164), 'yellow');
% plot(x_data(1:41), y_data(165:205), 'magenta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' Minimum location')
disp(['X  = ' , num2str(min_x_hat), ' ', num2str(min_x_alpha_l) , '  +',num2str(min_x_alpha_u)]);
disp(['Y  = ' , num2str(min_y_hat), ' ', num2str(min_y_alpha_l) , '  +',num2str(min_y_alpha_u)]);
disp(['Confidence interval= ', num2str(100*(1-alpha)), '%']);
disp(['Bootstrap samples= ', num2str(nboot)]);
disp(['Spline smoothing factor= ', num2str(smoothing_factor)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save x-position
save('min_x_star','min_x_star');
    
end

 end % regression_bootstrap function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_star = bootfun1(x, y)
global smoothing_factor;

% calculate f_star
f_star = csaps(x, y, smoothing_factor);

end % end bootfun1 function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_star = bootfun2(residuals)
global f_hat;
global x_data;
global y_data;
global smoothing_factor;


% calculate new y_star
y_star = residuals + y_data;

f_star = csaps(x_data, y_star, smoothing_factor);

end % bootfun2 function
    
