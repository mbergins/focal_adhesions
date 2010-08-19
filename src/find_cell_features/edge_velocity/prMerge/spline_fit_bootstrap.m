function [max_x_residuals] = spline_fit_bootstrap(x_data, y_data);
% spline_fit_bootstrap fits smoothing spline into data (x,y) with conf. int 
%
%
%       For details on bootstrap theory, standart error
%       and confidence intervalls see :
%       An Introduction to the Bootstrap
%       Robert J. Tibshirani  1993    
%                 
%
%
% SYNOPSIS      spline_fit_bootstrap(x_data, y_data)
%
% INPUT                :    x_data
%                           y_data
%      
% 
% OUTPUT               : 
%                           
%
% Matthias Machacek 10/20/06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear global f_hat;
clear global smoothing_factor;

global f_hat;
global x_data;
global y_data;
global smoothing_factor;

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

% spline smoothing factor 
% GTPases setting : smoothing_factor = 0.0001;
smoothing_factor = 0.8;
% smoothing_factor = 0.5;


xx_min = min(x_data);
xx_max = max(x_data);
xx_step = (xx_max - sign(xx_min)*abs(xx_min))/100;
if xx_step > 1
    xx_step = 1;
end
xx = xx_min:xx_step:xx_max;
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

% plot the smoothing spline with data points
figure
plot(xx, ppval(f_hat, xx));
hold on
plot(x_data, y_data, 'x');
%plot(x_data, y_hat_loess, ':');
%plot(x_data, y_hat_sgolay, ':');
title('Data with spline fit')'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now calculate the residuals using the estimation
residuals = y_data - ppval(f_hat, x_data);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%         Analyse residuals         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Get rms
rms = norm(residuals)/sqrt(length(residuals));
disp(['RMS =  ', num2str(rms)]);

figure
plot(x_data, residuals, 'x');
title('Residuals');
figure
hist(residuals,60);
title('Histogram of residuals');
% Estimate the density function of the residuals
[res_density, xi] = ksdensity(residuals);
figure
plot(xi, res_density);
title('Estimation of probability density fct from residuals');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     Compute all bootstrap regression estimates   %%%%%%%%%
% I commented the following two lines as it was plottin every booted sample
% and therefore crowding the fig. 5. Instead I added these tow lines after
% the following end.
figure  
hold on 
for i=1:nboot
    y_star(i,:) = ppval(bootstat(i), xx);
    % Extract position and value of peak for each f_star
    [max_y_star(i), max_index] = max(y_star(i,:));
    max_x_star(i) = xx(max_index);
    plot(xx,y_star(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%       Compare with estimate of f: f_hat
plot(xx, y_hat, 'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%       Analyse maxima    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get peak of f_hat
[max_y_hat, max_index] = max(y_hat);
max_x_hat = xx(max_index);

% compute residuals of the x and y positions of the peak
max_x_residuals = max_x_hat - max_x_star;
max_y_residuals = max_y_hat - max_y_star;

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
%%%%%%%%%%%%%%%%%%%%%     End analyse maxima     %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the individual data but colored
plot(x_data(1:41), y_data(1:41), 'black');
plot(x_data(1:41), y_data(42:82), 'blue');
plot(x_data(1:41), y_data(83:123), 'green');
plot(x_data(1:41), y_data(124:164), 'yellow');
plot(x_data(1:41), y_data(165:205), 'magenta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Plot results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot standart error
plot(xx, y_hat + se_star, 'c--');
plot(xx, y_hat - se_star, 'c--');

% plot confidence interval
plot(xx, y_hat + y_alpha_u, 'r--');
plot(xx, y_hat + y_alpha_l, 'r--');

plot(max_x_hat, max_y_hat,'black');

% plot standart error
myErrorbar(max_x_hat, max_y_hat,[se_max_x_star se_max_y_star], [se_max_x_star  se_max_y_star]);

% plot confidence interval
myErrorbar(max_x_hat, max_y_hat,[max_x_alpha_l max_y_alpha_u], [max_x_alpha_l  max_y_alpha_u]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' Peak location')
disp(['X  = ' , num2str(max_x_hat), ' ', num2str(max_x_alpha_l) , '  +',num2str(max_x_alpha_u)]);
disp(['Y  = ' , num2str(max_y_hat), ' ', num2str(max_y_alpha_l) , '  +',num2str(max_y_alpha_u)]);
disp(['Confidence interval= ', num2str(100*(1-alpha)), '%']);
disp(['Bootstrap samples= ', num2str(nboot)]);
disp(['Spline smoothing factor= ', num2str(smoothing_factor)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save x-position
save('max_x_star','max_x_star');


end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f_star = bootfun1(x, y)
global smoothing_factor;

% calculate f_star
f_star = csaps(x, y, smoothing_factor);

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_star = bootfun2(residuals)
global f_hat;
global x_data;
global y_data;
global smoothing_factor;

% calculate new y_star
y_star = residuals + y_data;

f_star = csaps(x_data, y_star, smoothing_factor);

end % function
    
