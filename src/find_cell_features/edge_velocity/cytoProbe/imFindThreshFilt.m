function [ans, thresh_val, max_bg_val, max_obj_val]=imFindThreshFilt(img_in, DEPTH, CONTR, varargin)
% IMFINDTHRESHFILT finds threshold separating image into background&object
% 
%             This routine can be used for images with a uniform background
%             and objects of one single intensity distribution. 
%             
%             Default settings are:
%             F_WINDOW=20;
%             F_SIGMA=0.6;
%             LOWER_CUT=0.005;
%             UPPER_CUT=0.99;
%             REL_RELEVANCE=0.6;
%             TOT_RELEVANCE=0.005;  
%
%             If the function fails ans = -1 is returned
%             otherwise ans = 1
%
% SYNOPSIS    [ans, thresh, max_bg, max_obj]=imFindThresh(img_in, DEPTH, CONTR)
%
% INPUT       img       : the image
%             DEPTH     : depth of the image
%             CONTR     : flag for control image display
% 
% OUTPUT      ans           : succes flag
%             thresh_val    : the threshold value
%             max_bg_val   : the mean value of the background
%             max_obj_val  : the mean value of the object
%                           
% DEPENDENCES   imFindThresh uses { extr_relevant_max
%                                 } 
%               imFindThresh is used by { imFindCellEdge
%                                 } 
%
% Matthias Machacek 25/09/03

%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if DEPTH == 8 | DEPTH == 10 | DEPTH == 12 | DEPTH == 14 | DEPTH == 16 | DEPTH == 255 ...
      | DEPTH == 1023 | DEPTH == 4095 | DEPTH == 16383 | DEPTH == 65535 
else 
   error('unsuported image depth, only images of 8,10,12,14,16 bit are accepted.');
end

%convert
switch DEPTH
   case 8
      DEPTH=2^8-1;
   case 10
      DEPTH=2^10-1;
   case 12
      DEPTH=2^12-1;
    case 14
      DEPTH=2^14-1;  
    case 16
      DEPTH=2^16-1;
end

l=length(varargin);
for i=1:2:l
    in_found=0;
    if strcmp(varargin{i},'f_window')
        F_WINDOW=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'f_sigma')
        F_SIGMA=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'lower_cut')
        LOWER_CUT=varargin{i+1}; 
        in_found=1; 
    elseif strcmp(varargin(i),'upper_cut')
        UPPER_CUT=varargin{i+1};  
        in_found=1; 
    elseif strcmp(varargin(i),'tot_relevance')
        REL_RELEVANCE=varargin{i+1};  
        in_found=1; 
    elseif strcmp(varargin(i),'rel_relevance')
        TOT_RELEVANCE=varargin{i+1};     
        in_found=1; 
    end
   
    if in_found == 0
        error_string = char(varargin(i));
        error(['Unknown input:   ' , error_string]);
    end        
end

% default filter parameter
if ~exist('F_WINDOW','var')
   F_WINDOW=10;
end
if ~exist('F_SIGMA','var')
   F_SIGMA=6;
end

%histogram acceptance levels in percentage of total pixel counts
if ~exist('LOWER_CUT','var')
   LOWER_CUT=0.005;
end
if ~exist('UPPER_CUT','var')
   UPPER_CUT=0.99;
end
if ~exist('REL_RELEVANCE','var')
   REL_RELEVANCE=0.05;%0.4;
end
if ~exist('TOT_RELEVANCE','var')
   TOT_RELEVANCE=0.005;
end
%%%%%%%%%%%%%%%%%%% End parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do a histogram binning
%put code in here

%calculate the histogram
[img_hist, hist_val] = imhist(img_in,DEPTH);


%take xx% percent of the pixels
total_pix=sum(img_hist);
%find 3% 
i=1;
while sum(img_hist(1:i)) < LOWER_CUT*total_pix
   i=i+1;
end
i_lower=i;
%find 97%
i=1;
while sum(img_hist(1:i)) < UPPER_CUT*total_pix
   i=i+1;
end
i_upper=i;

img_hist(1:i_lower)=0;
hist_val(1:i_lower)=0;
hist_length = length(img_hist);
img_hist(i_upper:hist_length)=[];
hist_val(i_upper:hist_length)=[];

% filter the histogram to smooth it
if DEPTH == 65535
    F_WINDOW=17;
    F_SIGMA=3.8;
elseif  DEPTH == 16383
    F_WINDOW=15;
    F_SIGMA=8;  
elseif  DEPTH == 255
    F_WINDOW=15;
    F_SIGMA=15;        
elseif  DEPTH < 16383
    F_WINDOW=15;
    F_SIGMA=8;    
end

f_sigma_fourier=1/(2*pi/F_SIGMA)/DEPTH;

%Data must have length more than 3 times filter order
%test if this is the case
if 3 * F_WINDOW >= length(img_hist)
    F_WINDOW = floor(length(img_hist) / 3);
end
w=gausswin(F_WINDOW,F_SIGMA);
img_hist_ff=filtfilt(w./sum(w),1,img_hist); 


if CONTR
    fig_hist_h = figure;
    subplot(2,1,1);plot(hist_val, img_hist);
    hold on
    subplot(2,1,1);plot(hist_val, img_hist_ff,'r');
    hist_fft=fft(img_hist);
    hist_spectrum = hist_fft.* conj(hist_fft) / length(hist_fft);
    frequ = linspace(0,0.5,floor(length(img_hist)/2));
    subplot(2,1,2);plot(frequ, hist_spectrum(1:floor(length(img_hist)/2)));
    h_axes = gca; 
    ylimits = get(h_axes(1),'YLim');
    line([f_sigma_fourier f_sigma_fourier], [0 ylimits(2)]);
    title('Image intensity power spectrum');   
end

% find local maximas
hist_ff_max=imregionalmax(img_hist_ff);
% find local minimas
hist_ff_min=imregionalmin(img_hist_ff);

%test maximas ans minimas on significance
[i_max j_max val_max]=find(hist_ff_max.*img_hist_ff);
[i_min j_min val_min]=find(hist_ff_min.*img_hist_ff);
    
%find the relevant maximas from all local maximas
[ind_max, ind_min]=extr_relevant_max(i_max, val_max, i_min, val_min, 'rel_relevance', REL_RELEVANCE, 'tot_relevance', TOT_RELEVANCE);

if ind_max == -99
    %no solution was found
    ans = -1;
    thresh_val  = 0;
    max_bg_val  = 0;
    max_obj_val = 0;
    return
end

if length(ind_max)<2
    disp('Only one object found!');
    ans=-1;
    hold off
elseif length(ind_min)==0
    disp('No minima (thus threshold found)');
    ans=-1;
    hold off
else
    ans=1;
    
    thresh = ind_min(1);
    max_bg = ind_max(1);
    max_obj= ind_max(2);

    %convert the intesity values from a [1 DEPTH] range
    %back to a normalized range [0 1]
    thresh_norm = ind_min(1)/DEPTH;
    max_bg_norm = ind_max(1)/DEPTH;
    max_obj_norm= ind_max(2)/DEPTH;
    
    thresh_val = ind_min(1)/DEPTH;
    max_bg_val = ind_max(1)/DEPTH;
    max_obj_val= ind_max(2)/DEPTH;  

    %calculate the variance of the gauss filter 
    %fourie transform
    f_sigma_fourier=1/(2*pi*F_SIGMA);

	if CONTR
        figure(fig_hist_h);
        subplot(2,1,1);plot(hist_val(ind_max'), img_hist_ff(ind_max'),'xr');
        subplot(2,1,1);plot(hist_val(ind_min'), img_hist_ff(ind_min'),'og');
        title('Image histogram');       
        hold off
	end
end



% %optional fitting of the intensity distribution by Gaussians
% if GAUSS_FIT
%     %eliminate zero elements from filtered histogram
%     mask=img_hist_ff>0;
%     [offset, j_temp, img_hist_ff_crop]= find(img_hist_ff);
%     img_hist_ff_crop=img_hist_ff_crop';
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % [ii jj]=find(img_hist);
%     % img_hist(ii(length(ii)):length(img_hist))=[];
%     % img_hist(1:ii(1))=[];
%     % t0=round(length(img_hist)/2);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %eliminate zero elements from original histogram
%     %img_hist_crop=img_hist.*mask;
%     %[offset, j_temp, img_hist_crop]= find(img_hist_crop);
%     %img_hist_crop=img_hist_crop';
%     img_hist_crop=img_hist(offset);
%     
% 	%variance
% 	x0(1)=1;
%     x0(2)=1;
% 	%mean
% 	x0(3)=ind_max(1)-offset(1);
%     x0(4)=ind_max(2)-offset(1);
% 	%scale
% 	x0(5)=img_hist_ff_crop(x0(3));
%     x0(6)=img_hist_ff_crop(x0(4));
% 	
% 	option=optimset('largescale','off','LevenbergMarquardt','on','display','none','showstatus','iterplus','gradobj','on','Diagnostics','off');
% 	[x,resnorm] = lsqnonlin(@GaussFunDiff,x0,[],[],option,img_hist_crop);
% 	
% 	[n m]=size(img_hist_crop);
% 	sigma1=x(1);
%     sigma2=x(2);
% 	m1=x(3);
%     m2=x(4);
% 	s1=x(5);
%     s2=x(6);
% 	
% 	i=1:m;
% 	g1 = s1*exp(-(i-m1).^2/(2*sigma1^2));
% 	g2 = s2*exp(-(i-m2).^2/(2*sigma2^2));
% 	if 1
% 	    figure
% 	    plot(g1+g2);
% 	    hold on
% 	    plot(g1);plot(g2);
% 	    plot(img_hist_ff_crop,'color','r');
%     end
% 	%get the minima of g1+g2
% 	hist_min= nonzeros(imregionalmin(g1+g2).*(g1+g2));
% 	hist_max= nonzeros(imregionalmax(g1+g2).*(g1+g2));
% 	[n_min_peaks, m]=size(hist_min);
% 	[n_max_peaks, m]=size(hist_max);
% 	
% 	%  1. find the minimum of g1+g2  and use it as initial value to find thresh   
% 	[i, j, val]=find(imregionalmin(g1+g2).*(g1+g2));
% 	%get the number of minimas
% 	nr_min=length(j);
% 	if nr_min == 3
% 	    thresh_init=j(2);
% 	elseif nr_min == 2
% 	    thresh_init=j(1);  
% 	else
% 	    error('unkown histogramm topology');
% 	end
% 	
% 	%  2. find the intersection of g1,g2
% 	x=0:0.01:length(g1);
% 	options = optimset('Display','none');  % Turn off Display
% 	thresh_cross = fsolve(@GaussFunMin,thresh_init,options,sigma1,sigma2,m1,m2,s1,s2);
%     
%     %convert it back
%     thresh=(thresh_cross+offset(1)+1)/DEPTH;
%     thresh_norm=thresh;
% 	
% 	%  3. find the thresh value based on a qual area criterion
% 	options = optimset('Display','none');  % Turn off Display
% 	thresh_area = fsolve(@AreaFun,thresh_init,options,sigma1,sigma2,m1,m2,s1,s2);
% 	
%     %convert it back
%     thresh=(thresh_area+offset(1)-1)/DEPTH;
%     thresh_norm=thresh;
%     
% 	%get the two mean values
% 	max_bg=find(imregionalmax(g1));
% 	max_obj=find(imregionalmax(g2));
% 
% 	figure
% 	plot(g1+g2);
% 	hold on
% 	plot(img_hist,'color','r'); 
% 	line([thresh_init;thresh_init],[0;10]);
% 	line([thresh_p;thresh_p],[0;10]);
% 	line([thresh_area_p;thresh_area_p],[0;10]);
% 	text(thresh_init,13,['thresh init= ',num2str(thresh_init)]);
% 	text(thresh_p     ,12,['thresh= ',num2str(thresh)]);
% 	text(thresh_area_p,10,['thresh area= ',num2str(thresh_area)]);
% 	title(['Histogramm, Filtered: ',num2str(GAUSS_F)]);
% end
%  
% function F = GaussFunDiff(x,data_points)
% % GAUSSFUNDIFF calc. diff. between data points and two gaussian
% % 
% %               
% % SYNOPSIS      F = GaussFunDiff(x,data_points)
% %
% % INPUT         x           : the 6 parameters of the two gaussian distributions
% %               data_points : the data_points to be approximated by the two
% %                               Gaussian
% % 
% % OUTPUT        F           : array with the difference between the data_points 
% %                             and the two gaussians
% %               GaussFunDiff is used by { imFindThresh }
% %
% % Matthias Machacek 25/09/03
% 
% sigma1=x(1);
% sigma2=x(2);
% m1=x(3);
% m2=x(4);
% s1=x(5);
% s2=x(6);
% 
% [n m]=size(data_points);
% 
% i=1:m;
% g1 = s1*exp(-(i-m1).^2/(2*sigma1^2));
% g2 = s2*exp(-(i-m2).^2/(2*sigma2^2));
% 
% F=data_points-g1-g2;
% 
% 
% 
% function F = AreaFun(x,sigma1,sigma2,m1,m2,s1,s2)
% % AREAFUN calc. diff. between data points and two gaussian
% % 
% %               
% % SYNOPSIS      F = AreaFun(x,data_points)
% %
% % INPUT         x           : the 6 parameters of the two gaussian distributions
% %               data_points : the data_points to be approximated by the two
% %                               Gaussian
% % 
% % OUTPUT        F           : array with the difference between the data_points 
% %                             and the two gaussians
% %               GaussFunDiff is used by { imFindThresh }
% %
% % Matthias Machacek 20/11/03
% 
% a1 = quad(@GaussEFun,    x, 1000,  [], [], sigma1,m1,s1);
% a2 = quad(@GaussEFun,    0,    x,  [], [], sigma2,m2,s2);
% 
% F=a1-a2;
% 
% function f = GaussEFun(x, sigma,m,s)
% f = s*exp(-(x-m).^2/(2*sigma^2));





