function [nx_n, ny_n]=prGetProtRegion(img, sp_x_n, sp_y_n, r_n, varargin)
% PRGETPROTREGION calculates the unit normal relative to the edge and a band mask
%        
%               This function calculates the unit normals along a edge
%               given by the two component spline at the paremeters given
%               by "r_n". The normal points away from the object assuming
%               that the object has a higher average intensity.
%               Additionally a band mask of depth PROT_DEPTH is calculated.
%
% SYNOPSIS    [nx_n, ny_n, band_mask]=prGetProtRegion(img, sp_x_n, sp_y_n, r_n)
%
% INPUT       img       : normalized grayscale image
%             sp_x_n    : edge x B-spline 
%             sp_x_n    : edge y B-spline 
%             r_n       : spline parameter
%
% 
% OUTPUT      nx_n      : x-component of the edge unit normal
%             ny_n      : y-component of the edge unit normal
%             band_mask : BW mask along the edge with a depth PROT_DEPTH
%   
%                           
% DEPENDENCES   prGetDisplacement uses { 
%                                       }
%
%               prGetDisplacement is used by { imEdgeTracker                                     
%                                       }
% Matthias Machacek 11/11/03

%  _nn variable at timestep t-1
%  _n  variable at timestep t


l=length(varargin);
for i=1:l
    if strcmp(varargin(i),'prot_depth')
        PROT_DEPTH=varargin{i+1};
    end
end

%default constants
if ~exist('PROT_DEPTH','var')
    PROT_DEPTH=30;
end

[n_img, m_img]=size(img);

knots_nr=length(r_n);

x_n=fnval(sp_x_n,r_n);
y_n=fnval(sp_y_n,r_n);

%derivative spline 
sp_dx_n = fnder(sp_x_n);
sp_dy_n = fnder(sp_y_n);
%derivatives at discrite locations
dx_n=fnval(sp_dx_n,r_n);
dy_n=fnval(sp_dy_n,r_n);
%normalize
l=sqrt(dx_n.^2+dy_n.^2);
dx_nn=dx_n./l;
dy_nn=dy_n./l;
%the normal unit vector (not oriented!!)
nx_n= dy_n;
ny_n=-dx_n;

%determine the object side of the edge
p1_x=round(x_n+2*nx_n);
p1_y=round(y_n+2*ny_n);

p2_x=round(x_n-2*nx_n);
p2_y=round(y_n-2*ny_n);

in_out1=0;
in_out2=0;

for i=1:length(p1_x)
    logic_exp1 = p1_x(i) >= 1     & p1_y(i) >= 1     & p2_x(i) >= 1     & p2_y(i) >= 1;
    logic_exp2 = p1_x(i) <= m_img & p1_y(i) <= n_img & p2_x(i) <= m_img & p2_y(i) <= n_img;
    if logic_exp1 & logic_exp2
        in_out1 = in_out1 + img(p1_y(i), p1_x(i));
        in_out2 = in_out2 + img(p2_y(i), p2_x(i));
    end
end

%assume that the background has a lower intensity
%the normal is pointin away from the object!
if in_out1 > in_out2
   nx_n= -nx_n;
   ny_n= -ny_n; 
end

%average normals
max_spline_par=r_n(length(r_n));
[nx_n_av, ny_n_av, x_n_av, y_n_av, pos] = prGetAvEdge(sp_x_n, sp_y_n,  max_spline_par, r_n, nx_n, ny_n, 'nr_sect', 20);

%re-nornamize the averaged normals
l_n_av = sqrt(nx_n_av.^2 + ny_n_av.^2);
nx_n_av = nx_n_av ./ l_n_av;
ny_n_av = ny_n_av ./ l_n_av;

c=x_n;
r=y_n;
j=length(nx_n_av)+1;
for i=1:length(nx_n_av)
    c(length(x_n)+i)= x_n_av(j-i) - PROT_DEPTH *nx_n_av(j-i);
    if c(length(x_n)+i) < 1 
        c(length(x_n)+i) = 1;
    elseif c(length(x_n)+i) > m_img
        c(length(x_n)+i) = m_img;
    end
    r(length(y_n)+i)= y_n_av(j-i) - PROT_DEPTH *ny_n_av(j-i);
    if r(length(x_n)+i) < 1 
        r(length(x_n)+i) = 1;
    elseif r(length(x_n)+i) > n_img
        r(length(x_n)+i) = n_img;
    end
end

%extract image segment
%band_mask=roipoly(img,c,r);

