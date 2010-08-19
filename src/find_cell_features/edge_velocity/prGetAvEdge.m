function [x_av, y_av, x_n_av, y_n_av, pos] = prGetAvEdge(sp_x, sp_y, max_spline_par, p_vec, x_vec, y_vec, varargin)
% PRGETAVEDGE averages vectors assiciated along a spline
% 
%             Given a spline and values along that spline, this function
%             devides the spline in a number of segments and calculates the
%             average of this vectors in each segment.
%             It is assumed that the spline parameter starts with 1 and
%             ends with max_spline_par. If there is no mean value found for
%             a segment the values are set to NaN.
%             
%
% SYNOPSIS    [x_av, y_av, x_n_av, y_n_av, pos] = prGetAvEdge(sp_x, sp_y, max_spline_par, p_vec, x_vec, y_vec, varargin)
%
% INPUT       sp_x              : x-spline
%             sp_y              : y-spline
%             max_spline_par    : max value of the spline parameter
%             p_vec             : parameter values of the vectors
%             x_vec             : x-comp. of the vectors
%             y_vec             : y-comp. of the vectors
% 
% OUTPUT      x_av              : x-comp of averaged vector
%             y_av              : y-comp of averaged vector
%             x_n_av            : x-comp of averaged vector position
%             y_n_av            : x-comp of averaged vector position
%             pos               : parameters of the averaged vector position
%                           
% DEPENDENCES   prGetAvEdge uses { 
%                                 } 
%               prGetAvEdge is used by { imEdgeTracker
%                                 } 
%
% Matthias Machacek 25/11/03

%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);

for i=1:l
   if strcmp(varargin(i),'nr_sect')
      NR_SECT=varargin{i+1};    
   end
end

if ~exist('NR_SECT','var')
   NR_SECT=10;
end
%%%%%%%%%%%%%%%%%%% End parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the spline parameter interval based
% on given number of sections
interval = (max_spline_par-1)/NR_SECT;

% create a vector of section boundaries
pos = 1 : interval : max_spline_par;

% average vectors along the edge
i=1;
for j=1:NR_SECT
    nr=0;
    x_av(j)=0;
    y_av(j)=0;
    while i < length(p_vec) && p_vec(i) < pos(j+1)   %before length(p_vec)!!!
        x_av(j) =  x_av(j) + x_vec(i);
        y_av(j) =  y_av(j) + y_vec(i); 
        nr=nr+1;
        i=i+1;
    end
    if nr > 0 
        % averaged values
        x_av(j)=x_av(j)/nr;
        y_av(j)=y_av(j)/nr;     
    else
        x_av(j)= NaN;
        y_av(j)= NaN;  
    end
    
    % spline parameters of the averaged vector position
    
    % positions of the averaged values on the spline
    x_n_av(j)=fnval(sp_x,(pos(j+1)-pos(j))/2+pos(j));
    y_n_av(j)=fnval(sp_y,(pos(j+1)-pos(j))/2+pos(j));
    % quiver(x_n_av, y_n_av, x_normal_out_av(j), y_normal_out_av(j),10,'r');
end