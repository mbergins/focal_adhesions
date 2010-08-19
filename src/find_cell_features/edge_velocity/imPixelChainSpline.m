function [sp_x_s, sp_y_s]=imPixelChainSpline(pixel_list, varargin)
% IMPIXELCHAINSPLINE given one conected set pixels it calc. two-parametric spline
% 
%              The default settings are:
%              TOLERANCE=20;
%              SMOOTH=1;
%              NORMALIZE=0;
%
%
% SYNOPSIS    [sp_x_s, sp_y_s]=imPixelChainSpline(pixel_list)
%
% INPUT       pixel list   : SORTED set of pixels
%             
% 
% OUTPUT      sp_x_s       : the spline for the x-coordinate
%             sp_y_s       : the spline for the y-coordinate
%                           
% DEPENDENCES   imPixelChainSpline is used by { imEdgeTracker
%                                       }
%
% Matthias Machacek 10/09/03


%%%%%%%%%%%%%%   Parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);

for i=1:l
   if strcmp(varargin{i},'tolerance')
      TOLERANCE=varargin{i+1};
   end
   if strcmp(varargin{i},'smooth')
      SMOOTH=1;
   end   
end

if ~exist('TOLERANCE','var')
   TOLERANCE=40;
end
if ~exist('SMOOTH','var')
   SMOOTH=1;
end
if ~exist('NORMALIZE','var')
   NORMALIZE=0;
end
%%%%%%%%%%%%%%   End parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[edge_l, dim]=size(pixel_list);
if NORMALIZE
   dl=1/(edge_l-1);
   i=0:dl:1;
   i=i';
else
   i=1:edge_l;
   i=i';
end

if SMOOTH
   %smoothing spline
   sp_x_s = spaps(i,pixel_list(:,1),TOLERANCE);
   sp_y_s = spaps(i,pixel_list(:,2),TOLERANCE);
else
	%least square spline
	%number of spline segments
	seg_nr=floor(edge_l/10);
	sp_x = spap2(seg_nr,4,i,pixel_list(:,1));
	sp_y = spap2(seg_nr,4,i,pixel_list(:,2));
end


