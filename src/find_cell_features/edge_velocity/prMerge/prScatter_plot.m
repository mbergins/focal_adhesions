function prScatter_plot(edge_parameters, merg_parameters, post_parameters, var1, var2, var1_name, var2_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyze the scatter plots                                     %
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
mkdir([RESULT_DIR filesep 'scatter' filesep]);
RESULT_DIR = [RESULT_DIR filesep 'scatter' filesep];

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


% Shift variable 2 by the specified number of time points
if SCO_TIME_SHIFT > 0
    var2 = circshift(var2,[0 SCO_TIME_SHIFT]);
    var2(:,1:SCO_TIME_SHIFT)=[];
    var1(:,1:SCO_TIME_SHIFT)=[];
elseif SCO_TIME_SHIFT < 0
    var2 = circshift(var2,[0 SCO_TIME_SHIFT]);
    var2(:,end+SCO_TIME_SHIFT:end)=[];
    var1(:,end+SCO_TIME_SHIFT:end)=[];
end
    
% Put matrix into array
var1_array = reshape(var1,1,[]);
var2_array = reshape(var2,1,[]);

% remove zeros
index_z = find(~var2_array);
var1_array(index_z) = [];
var2_array(index_z) = [];

% do principal component analysis on normalized variables
stdr_var1 = std(var1_array);
stdr_var2 = std(var2_array);
[V,D] = eigs(cov(var1_array./stdr_var1, var2_array./stdr_var2));
mean_var1 = mean(var1_array);
mean_var2 = mean(var2_array);

h_scatter1 = figure;
%plot(var1_array ,var2_array,'.');
bins = 25; 
CMAP = flipud(gray(256));
%CMAP = flipud(jet(256));
scattercloud(var1_array ,var2_array, bins, 1, 'k.', CMAP);
xlabel(var1_name);
ylabel(var2_name);
title('Original variables');


h_scatter2 = figure;
%plot(var1_array./stdr_var1 ,var2_array./stdr_var2,'.');
bins = 25; 
CMAP = flipud(gray(256));
%CMAP = flipud(jet(256));
scattercloud(var1_array./stdr_var1 ,var2_array./stdr_var2, bins, 1, 'k.', CMAP);
% plot principal components
hold on
plot(D(1,1).*[0,V(1,1)]+mean_var1./stdr_var1,   D(1,1).*[0, V(2,1)]+mean_var2./stdr_var2,'r-');
plot(D(2,2).*[0,V(1,2)]+mean_var1./stdr_var1,   D(2,2).*[0, V(2,2)]+mean_var2./stdr_var2,'g-');
xlabel(var1_name);
ylabel(var2_name);
title(['Variables are normalized ! '' Ratio of the 1st to 2nd principal component: ' num2str(D(1,1)/D(2,2))]);

result_file_name = ['scatter_', var1_name, var2_name]
%hgsave(h_scatter,[RESULT_DIR result_file_name]);
%print(h_scatter, [RESULT_DIR result_file_name],'-depsc2','-tiff');
%print(h_scatter, [RESULT_DIR result_file_name],'-dtiff');

return

% % Robust linear fitting
% [b_var1_var2, statsb_var1_var2] = robustfit(var1_array, var2_array);
%     
% %[p,S,mu] = polyfit(pr,sc,1);
% %[y,delta] = polyconf(p,pr,S);
% %sc_poly_fit = polyval(p,pr);
% % try as-hoc method
% var1_min = min(var1);
% var1_max = max(var1);
% var1_step = (max(var1)-min(var1))/8;
% ii=1;
% for i = var1_min : var1_step : pr_max
%     [i1, i2, var1_val] = find(var1 > i & var1 < (i + var1_step));
%     av_var2(ii)  = mean(var2(i2));
%     std_var2(ii) = std(var2(i2));
%     ii=ii+1;
% end



% hold on
% plot(pr_min : pr_step : pr_max, av_re,'r');
% plot(pr_min : pr_step : pr_max, av_re+std_re,'r--');
% plot(pr_min : pr_step : pr_max, av_re-std_re,'r--');
% plot(pr_depoly,re_depoly,'g+');
% plot([pr_min pr_max],[b_pr_re(1) + b_pr_re(2)*pr_min, b_pr_re(1) + b_pr_re(2)*pr_max],'r');
% 
% plot(av_pr_poly, av_re_poly, 'bo');
% plot(av_pr_depoly, av_re_depoly, 'bo');



    
%     poly_mask = find(sc > 0);
%     depoly_mask = find(sc < 0);
%     pr_poly     = pr(poly_mask);
%     pr_depoly   = pr(depoly_mask);
%     re_poly     = re(poly_mask);
%     re_depoly   = re(depoly_mask); `1 
%     
%     % get center of poly and depoly
%     av_re_poly = mean(re_poly,2);
%     av_pr_poly = mean(pr_poly,2);
%     av_re_depoly = mean(re_depoly,2);
%     av_pr_depoly = mean(pr_depoly,2);
    

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%               do cluster analysis      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    verbose = 1;
    %accuracy
    th = 10^(-4);
    regularize = 0;
    covoption = 0;
    s_vec = [pr',re', sc']';
    [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(s_vec, 1, 10,...
        regularize, th, covoption, [], [], verbose);

    %     figure(h_var_activity_space);
    %     hold on
    %     for comp=1:length(bestpp)
    %         elipsnorm(bestmu([1,2],comp),bestcov([1,2],[1,2],comp),2)
    %     end
    %     drawnow

    for k = 1:bestk
        rho(k) = bestcov(1,2,k)/(sqrt(bestcov(1,1,k))*sqrt(bestcov(2,2,k)));
        a(k) = atan(2* bestcov(1,2,k)/(bestcov(1,1,k)-bestcov(2,2,k)))/2;

        p1(k)=sqrt(  bestcov(1,1,k) * bestcov(2,2,k)*(1-rho(k)^2)/...
            (bestcov(2,2,k)*cos(a(k))^2 - ...
            2* rho(k)*sqrt(bestcov(1,1,k))*sqrt(bestcov(2,2,k))*sin(a(k))*cos(a(k)) +...
            bestcov(1,1,k)*sin(a(k))^2));
        p2(k)=sqrt(  bestcov(1,1,k) * bestcov(2,2,k)*(1-rho(k)^2)/...
            (bestcov(2,2,k)*sin(a(k))^2 - ...
            2* rho(k)*sqrt(bestcov(1,1,k))*sqrt(bestcov(2,2,k))*sin(a(k))*cos(a(k)) +...
            bestcov(1,1,k)*cos(a(k))^2));
    end


    %calculate the probability of each data point
    for k = 1:bestk
        covar = bestcov(:,:,k);
        m = bestmu(:,k);
        y(:,k) = multinorm(s_vec,m,covar);

        %ff = ((2*pi*(var+realmin))^(-1/2));
        %y(:,k) = bestpp(k) * ff * exp((-1/(2*var))*(s_vec-m).^2);
    end

    %find the highest probability for each data point
    [val, cluster_index] = max(y,[],2);

    %sort the index according to the intensity
    [cluster_c  cluster_c_index]= sort(bestmu', 1);
    cluster_index_sort = zeros(length(cluster_index),1);
    for i=1:bestk
        cluster_index_sort = cluster_index_sort + (i * (cluster_index == cluster_c_index(i)));
    end

    % check for zeros in s_vec and a small irrelevant value
    z1 = s_vec(1,:)==0;
    s_vec(1,:)=s_vec(1,:)+z1*0.00001;
    z2 = s_vec(2,:)==0;
    s_vec(2,:)=s_vec(2,:)+z2*0.00001;

    for i=1:bestk
        [d1,d2,cluster{i,1}(1,:)] = find((cluster_index == i)'  .* s_vec(1,:));
        [d1,d2,cluster{i,1}(2,:)] = find((cluster_index == i)'  .* s_vec(2,:));
    end

    col=['r','g','b','c','y','r'];
    h_cluster = figure;
    hold on
    for i=1:bestk
        plot(cluster{i,1}(1,:), cluster{i,1}(2,:),'.','MarkerEdgeColor',col(i));
    end
    for i=1:bestk
        plot(bestmu(1,i),bestmu(2,i),'s','MarkerEdgeColor','b');
        line([bestmu(1,i)-1.5*sqrt(bestcov(1,1,i)) ,bestmu(1,i)+1.5*sqrt(bestcov(1,1,i)) ],...
            [bestmu(2,i)-1.5*sqrt(bestcov(2,2,i)) ,bestmu(2,i)+1.5*sqrt(bestcov(2,2,i)) ],'Color','k');
        %             line([bestmu(1,i)-p1(i) ,bestmu(1,i)+p1(i) ],...
        %                  [bestmu(2,i)-p2(i) ,bestmu(2,i)+p2(i)],'Color','k');
        text(bestmu(1,i),bestmu(2,i),['   (',num2str(bestmu(1,i)),',  ', num2str(bestmu(2,i)),')']);
        %elipsnorm(bestmu([1,2],i),bestcov([1,2],[1,2],i),2)
    end
    xlabel('Protrusion [nm/s]');
    ylabel('Activity');
    hgsave(h_cluster,[RESULT_DIR 'prot_activity_cluster.fig']);
    print(h_cluster, [RESULT_DIR 'prot_activity_cluster.eps'],'-depsc2','-tiff');
    print(h_cluster, [RESULT_DIR 'prot_activity_cluster.tif'],'-dtiff');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

