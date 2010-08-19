function [data_filtered, data_noise] = prFilterTimeSeries(data, FILTER_METHOD, windowSize, spline_s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyse the time series                                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DO_ANALYSIS = 0;
if DO_ANALYSIS
    significance = 0.005;
    
    % Checking for linear trend
    [HlinearTrend,pValueLinTrend,errFlag] = differenceSignTest(traj,significance);

    % Check if series is IID
    [Hidd,pValue,errFlag] = turningPointTest(traj,significance);
    
    % test if time series is IID by looking at its autocorrelation function.
    [HiidAC,pValue,errFlag] = portmanteau(traj,maxLag,significance)


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data_mean = mean(data,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detrending by difference operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD == 2
    delta_t = 1;

    data_detrend = diff(data,delta_t,2);
    data_detrend(:,end+1) = data_detrend(:,end);
    %data_trend = data - circshift(data_detrend, [0 1]);
    data_trend = data - data_detrend;
    
    % take the average of the detrended data
    av_data_detrend = mean(data_detrend,1);
    av_data_trend = mean(data_trend,1);
    
    
    % variance (sigma^2)
    var_trend   = var(data_trend,1,2);
    var_detrend = var(data_detrend,1,2);
    
    data_mean_detrend = diff(data_mean,delta_t,2);
    data_mean_detrend(:,end+1) = data_mean_detrend(:,end);
    mean_trend = data_mean - data_mean_detrend;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter by polynomial approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD == 3 | FILTER_METHOD == 4
    w1 = 1:size(data,2);
    w  = repmat(w1,size(data,1),1);
    
    if FILTER_METHOD == 3
        p = polyfit(w,data,4);
    elseif FILTER_METHOD == 4
        p = polyfit(w,data,5);
    end
    
    data_trend = polyval(p,w);
    data_detrend = data - data_trend;
    
    % take the average of the detrended data
    av_data_detrend = mean(data_detrend,1);
    av_data_trend = mean(data_trend,1);
    
    % variance (sigma^2)
    var_trend   = var(data_trend,1,2);
    var_detrend = var(data_detrend,1,2);
    
    data_mean_trend = mean(data_trend,1);
    data_mean_detrend = data_mean - data_mean_trend;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter by moving average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD == 5
    % initial condition
    %z = data(:,1:windowSize-1);
    %zi = flipdim(z,2)';
    z = data(:,2);
    zi = repmat(z,1,windowSize-1)'; 
    
    data_filtered = filter(ones(1,windowSize)/windowSize,1,data, zi, 2);
    data_noise = data - data_filtered;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detrending by Gauss filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD == 6
    
    % windowSize is the standart deviation sigma
    
    sig = windowSize;
    base = 3*windowSize;
    if base < 3
        base = 3;
    elseif ~isinteger(base) 
        base = ceil(base);   
    end

    filter_shape = fspecial('gaussian', [1 base], sig);
    %filter_shape = gausswin(sig,a);
    %data_filtered = imfilter(data, h,'symmetric');
    
    seg_nr = size(data,1);
    for i=1:seg_nr
        data_filtered(i,:) = filtfilt(filter_shape, 1, data(i,:));
        %data_filtered(i,:) = filter(filter_shape, 1, data(i,:));
    end
    data_noise = data - data_filtered;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detrending by Spline Interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD == 7
    
    tol = 20;
    
    % temporal fitting
    %x=1:size(data,2);
    
    % tempo-,spatial fitting
    x = {1:size(data,1),1:size(data,2)};
    
    %x=repmat(x,size(data,1),1);
    
    %sp = spaps(x,data,tol);
    
    %p=0 total smooth, p=1 interpolation
    %p = 0.4 .* ones(size(data,1),1);
    pp = csaps(x,data,spline_s); 
    
    data_filtered = fnval(pp,x);
    data_noise    = data - data_filtered;
end
