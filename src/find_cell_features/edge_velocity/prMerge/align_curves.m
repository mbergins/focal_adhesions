function [av_curve, dt_curves, av_dt_curves] = align_curves(curves, method)
% ALIGN_CURVE aligns a set of curves according to a certain criterion 
% 
%            
%
% SYNOPSIS        [aligned_curves, shifts, average_shift] = align_curves(curves, method)
%
% INPUT           curves:       nxm matrix where n is the number of curves
%                               and m the numer of elements in a curve
%                 method:       flag specifying the alignment criterion
% 
% OUTPUT          aligned_curves:
%                 shifts:
%                 average_shift: 
%                           
% DEPENDENCES     align_curves (                               
%                                   }
%                 align_curves{}
%
% Matthias Machacek 05/31/06


%    method = 0 : cross correlation based alignment
%    method = 1 : maximum value based alignment
%    method = 2 : minimum value based alignment


% plot correlation from all segments
% figure,waterfall(curves);
% xlabel('Time shift');
% ylabel('Curve');
% zlabel('Element');


n_curves   = size(curves,1);
n_elements = size(curves,2);     

% test curves for statistical significance

% Calculate the significance bounds
% 99% probability -> 2.58
% 95% probability -> 1.96
% 90% probability -> 1.645

upper_sig_bound = (n_elements)^(-0.5) * 1.96;
lower_sig_bound = -upper_sig_bound;

% test for significance of positive values
%hits = curves > upper_sig_bound;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% Step 1: get dt between curves:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%   Get shift with cross correlations  %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get mean curve as reference shape 
    ref_curve = mean(curves,1);
    lag_alignment = 5;
    for i_seg=1:n_curves
        curve_corr(i_seg,:) = xcov(ref_curve, curves(i_seg,:), lag_alignment,'coeff');
    end
    % Step 2: get the maximas from curve_shift to obtain the shifts
    for i_seg=1:n_curves
        [v,i] = max(curve_corr(i_seg,:));
        % check if it is a valid maxima
        if i>1 && i < lag_alignment;
            dt_curves(i_seg) = i;
        else
            dt_curves(i_seg) = NaN; 
        end
    end
    dt_curves = dt_curves-lag_alignment+1;
    figure,plot(dt_curves);
    xlabel('Sampling window number');
    ylabel('Shift correction');
    
elseif method == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%   Get shift with maximas             %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for i_seg=1:n_curves
        [v,i] = max(curves(i_seg,:));
        max_curves(i_seg) = i;
    end
    dt_ref = mean(max_curves);
    for i_seg=1:n_curves
        dt_curves(i_seg) = dt_ref - max_curves(i_seg);
    end
    figure,plot(dt_curves);
    xlabel('Sampling window');
    ylabel('Time shift');
elseif  method == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%   Get shift with maximas             %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
    
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3: align the curves                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average shift
av_dt_curves  = mean(dt_curves);

for i_seg=1:n_curves
    aligned_curves(i_seg,:) = circshift(curves(i_seg,:), [0, round(dt_curves(i_seg))]);
end
if 0
    figure,waterfall(aligned_curves);
    xlabel('Time shift');
    ylabel('Segment');
    zlabel('Correlation');
end

av_curve = mean(aligned_curves,1);
figure,plot(av_curve);
% compare
hold on
plot(mean(curves,1),'--');
legend('Aligned mean','Conventional mean');

% figure, plot(aligned_curves(5,:),'--');
% hold on
% plot(curves(5,:));
% plot(curves(1,:),'r');
% egend('Aligned curve nr 5','Conventional mean');


