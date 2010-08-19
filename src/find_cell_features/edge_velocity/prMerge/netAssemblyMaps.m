function [score_map, score_map_sign, num_scores] = netAssemblyMaps(SCORE, time, delta_t_av, imgSize, contr, varargin)
% NETASSEMBLYMAPS extracts, averages, filters and interpolates the scores
% 
%               SCORES has the following structure:
%               SCORES=[m,4], where
%               SCORES[m,1] is the time
%               SCORES[m,2] is the x-coordinate
%               SCORES[m,3] is the y-coordinate
%               SCORES[m,4] is the score              
%
% SYNOPSIS    netAssemblyMaps(SCORE, delta_t_av, imgSize, outputFileName)
%
% INPUT       SCORE             : the structure containing the score information
%             delta_t_av        : number of time steps used for averaging
%             imgSize           : image size
%             outputFileName    : outpu file name
% 
% OUTPUT      score_map         : map with the scores
%             num_scores        : number of measured scores
%                           
% DEPENDENCES   netAssemblyMaps uses {    
%                                       }
%
%               netAssemblyMaps is used by { prAlpha
%                                           }
%
% Matthias Machacek 12/04/03

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);
for i=1:l
    if strcmp(varargin(i),'filter_var')
        FILTER_VAR=varargin{i+1};
    elseif strcmp(varargin(i),'filter')
        FILTER=varargin{i+1}; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('FILTER_VAR','var')
    %filter width of the gauss
    FILTER_VAR = 5;
end
if ~exist('FILTER','var')
    % filter width of the gauss
    FILTER = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract scores
scores=SCORE(find(SCORE(:,1)>=time & SCORE(:,1)<=time+ delta_t_av -1),:);

% Initialize maps
score_map_sign = zeros(imgSize(1),imgSize(2),2);
score_map      = zeros(imgSize(1),imgSize(2));

% number of measurments
num_scores = size(scores,1);

% Put scores
for j=1:num_scores
    switch sign(scores(j,4))
        case 1
            score_map_sign(scores(j,2),scores(j,3),1)=score_map_sign(scores(j,2),scores(j,3),1) +    scores(j,4);
        case -1
            score_map_sign(scores(j,2),scores(j,3),2)=score_map_sign(scores(j,2),scores(j,3),2) +abs(scores(j,4));
        otherwise
            error('The score is zero.');
    end
    score_map(scores(j,2),scores(j,3)) = score_map(scores(j,2),scores(j,3))+scores(j,4);
end

if FILTER
    % Filter
    img_score_map(:,:,1)=Gauss2D(img_score_map(:,:,1),FILTER_VAR);
    img_score_map(:,:,2)=Gauss2D(img_score_map(:,:,2),FILTER_VAR);
    score_map(:,:)=Gauss2D(score_map(:,:),FILTER_VAR);
    if contr
        figure
        imshow(img_score_map,[]);
        title('The non filtered scores');
    end    
end

% i removed that becaue it does not make any sense. ?
%img_score_map=img_score_map./max(img_score_map(:));

% Display
if contr
    figure
    imshow(img_score_map,[]);
    title('Filtered scores');
end