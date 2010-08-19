function [seg_img, varargout]=imClusterSeg(img_in, CONTR, varargin)
% IMKMEANSSEG image segmentation with clustering (kmeans or EM)
% 
%             This function does a segmentation of the image using the 
%             kmeans or the EM algorithm. For the kmeans the number of
%             cluster (objects) has to be specified (k_cluster), for the 
%             EM algorithm the maximumal and minimal number of clusters
%             has to be specified (k_min, k_max) the function finds then
%             the optimal number of clusters. 
%             The EM algorithm can be slow, if so use image binning (2 to
%             4), and if running the function on image stacks use previous
%             results (pp, mu) as innitial results for the next image!
%             Watch out, if initial (pp, mu) are given their dimension must
%             match k_max!
%             The format of the segment image is as following:
%             I=[1..k_clusters] where I=1 corresponds to the class with the
%             lowes intensity and I=k_cluster to the class with the highest
%             intensity.
%             Binning can be used to speed up the function by binning with
%             a factor 0 to n. For the kmeans algorithm a coarsing occurs,
%             for the EM algorithm NO coarsing occurs!
%
%             If there is an error seg_img = -99 is returned
%  
%             Optional arguments are:
%               'method': 'em' or 'kmeans',
%               'binning': [0..n..N] where n is an integer
%             kmeans specific:
%               'mu0': initinal values, [] is also a valid input
%               'k_cluster': the number of clusters
%               'distance': 'sqEuclidean', 'cityblock' or 'cosine'
%             EM specific:
%               'k_min': minimum number of clusters
%               'k_max': maximum number of clusters
%               'p0': initinal values for individual cluster probability, [] is also valid
%               'mu0': initinal values for cluster mean values, [] is also a valid input
%
%
%
%
% SYNOPSIS    [seg_img, obj_val, varargout]=imClusterSeg(img_in, CONTR, varargin)
%
% INPUT       img_in    : the image
%             CONTR     : flag for control image display
% 
% OUTPUT      seg_img    : segmented image, intensities equal cluster #
%             optional   
%             p          : expectation of each cluster (only for EM)
%             mu         : the mean intensity of each cluster (sorted)
%             
%
%                           
% DEPENDENCES   imFindThresh uses { kmeans,
%                                   mixtures4
%                                 } 
%               imFindThresh is used by {
%                                 } 
%
% Matthias Machacek 02/10/04

%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);
for i=1:2:l
    in_found=0; 
    if strcmp(varargin{i},'method')
        METHOD=varargin{i+1};
        in_found=1;      
    elseif strcmp(varargin{i},'k_cluster')
        K_CLUSTER=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'k_min')
        K_MIN=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'k_max')
        K_MAX=varargin{i+1};  
        in_found=1; 
    elseif strcmp(varargin(i),'p0')
        P0=varargin{i+1};  
        in_found=1; 
    elseif strcmp(varargin(i),'mu0')
        MU0=varargin{i+1};      
        in_found=1; 
    elseif strcmp(varargin(i),'distance')
        DISTANCE=varargin{i+1};  
        in_found=1; 
    elseif strcmp(varargin(i),'binning')
        BINNING=varargin{i+1};       
        in_found=1; 
    elseif strcmp(varargin(i),'emptyaction')
        EMPTYACTION=varargin{i+1};  
        in_found=1; 
    end
    
    if in_found == 0
        error_string = char(varargin(i));
        error(['Unknown input:   ' , error_string]);
    end       
end

%%%%%%%%%  general parameter %%%%%%%%%%%%%%
if ~exist('METHOD','var')
    METHOD = 'em';
end
if ~exist('BINNING','var')
    BINNING = 0;
end

%%%%%%%%%%%  kmeans parameter %%%%%%%%%%%%%
if ~exist('K_CLUSTER','var')
   K_CLUSTER=2;
end
if ~exist('EMPTYACTION','var')
   EMPTYACTION='singleton';
   %'error'
   %'drop'
   %'singleton'
end
if ~exist('DISTANCE','var')
   DISTANCE='cityblock';
   %'cityblock'
   %'sqEuclidean'
   %'cosine'
end
%%%%%%%%%%   EM parameter %%%%%%%%%%%%%%%%
if ~exist('K_MIN','var')
   K_MIN=2;
end
if ~exist('K_MAX','var')
   K_MAX=4;
end
%the probability of the individual normal distributions
if ~exist('P0','var')
   P0=[];
end
%the mean values of the individual normal distributions
if ~exist('MU0','var')
   MU0=[];
end
%%%%%%%%%%%%%%%%%%% End parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_img_org m_img_org] = size(img_in);

%do image binning
if BINNING > 0
    img_bin = img_in(1:BINNING:n_img_org, 1:BINNING:m_img_org);
else 
    img_bin =  img_in;
end

%size of the image 
[n_img m_img] = size(img_bin);

%reshape image into vector
img_in_vec = reshape(img_bin, n_img*m_img, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  clustering with kmeans  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(METHOD,'kmeans')
    if ~isempty(MU0)
        if size(MU0) ~= K_CLUSTER
            h = warndlg('The size of your initial values for the clusters is different from your k_cluster','Warning!')    
        end
        %start clustering with initial values
        [cluster_index, cluster_c, cluster_d_sum]  =  kmeans(img_in_vec, K_CLUSTER, 'Distance', DISTANCE, 'start', MU0, 'emptyaction', EMPTYACTION); 
    else
        %start clustering without initial values
        [cluster_index, cluster_c, cluster_d_sum]  =  kmeans(img_in_vec, K_CLUSTER, 'Distance', DISTANCE, 'emptyaction', EMPTYACTION);
    end
    
    %sort the index according to the intensity
    [cluster_c  cluster_c_index]= sort(cluster_c, 1);
    cluster_index_sort = zeros(length(cluster_index),1);
    for i=1:K_CLUSTER
       cluster_index_sort = cluster_index_sort + (i * (cluster_index == cluster_c_index(i)));
    end
    
    %put vector back into matrix
    seg_img = reshape(cluster_index_sort, n_img, m_img);
    
    %resize image to original size
    if BINNING > 0
        seg_img = imresize(seg_img, [n_img_org, m_img_org], 'bicubic');   
%       for some reason this did not work
%         %determine the separating intensities
%         thresh_cluster_val(1) = 0;
%         for i=2:K_CLUSTER
%             thresh_cluster_val(i) = (cluster_c(i)+cluster_c(i-1))/2;
%         end
%         seg_img = zeros(n_img_org, m_img_org);
%         for i=1:K_CLUSTER-1
%             seg_img =  seg_img + i * ((img_in > thresh_cluster_val(i)) & (img_in < thresh_cluster_val(i+1)));
%         end
%         seg_img = seg_img + K_CLUSTER * (img_in > thresh_cluster_val(K_CLUSTER));
    end
    

    
    varargout{1} = [];
    varargout{2} = cluster_c;    
    
elseif strcmp(METHOD,'em')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  clustering with EM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %check first of the initial values have the right dimension
    if ~isempty(P0)
        if size(P0) ~= size(MU0)
            seg_img = -99;
            obj_val = -99;
            return 
        end
    end
    if ~isempty(MU0)
        if size(MU0) ~= size(MU0)
            seg_img = -99;
            obj_val = -99;
            return 
        end
    end    
    
    if CONTR
        verbose = 1;
    else
        verbose = 0;         
    end
    
    %accuracy
    th = 10^(-4);
    regularize = 0;
    covoption = 0;
    [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(img_in_vec', K_MIN, K_MAX,...
                                                         regularize, th, covoption, P0, MU0, verbose);
   
                                     
    %reshape image into vector
    img_in_vec = reshape(img_in, n_img_org*m_img_org, 1);                                                  
    %calculate the probability of each data point
    for k = 1:bestk
        var = bestcov(:,:,k);
        m = bestmu(:,k);
        ff = ((2*pi*(var+realmin))^(-1/2));
        y(:,k) = bestpp(k) * ff * exp((-1/(2*var))*(img_in_vec-m).^2);   
    end
    
    %find the highest probability for each data point
    [val, cluster_index] = max(y,[],2);     
    
    %sort the index according to the intensity
    [cluster_c  cluster_c_index]= sort(bestmu', 1);
    cluster_index_sort = zeros(length(cluster_index),1);
    for i=1:bestk
       cluster_index_sort = cluster_index_sort + (i * (cluster_index == cluster_c_index(i)));
    end
    
    %put vector into matrix shape
    seg_img = reshape(cluster_index_sort, n_img_org, m_img_org);

    %sort also the cluster expectation pp in the same way as the intensities 
    for j = 1:bestk
        bestpp_sort(j) = bestpp(cluster_c_index(j)); 
    end
    varargout{1} = bestpp_sort;
    varargout{2} = cluster_c';
end

if CONTR
    figure
    imshow(seg_img,[]);
end
