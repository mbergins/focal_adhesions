function [p_nn, disp, r_n, status]=prGetDisplacement(sp_x_nn, sp_y_nn, sp_x_n, sp_y_n, r_nn)
% PRGETDISPLACEMENT calculates displacement between two two-parametric (x(r),y(r)) splines
%                   here the displacement between the spline at t and the spline at 
%                   t+1 is defined as the distance measured at a 
%                   normal angle relative to the spline t
% 
% SYNOPSIS    [p_nn, disp]=prGetDisplacement(sp_x_nn, sp_x_nn, sp_x_n, sp_y_n, r_nn)
%
% INPUT       sp_x_n    : x-spline at time step t
%             sp_y_n    : y-spline at time step t
%             sp_x_nn   : x-spline at time step t-1
%             sp_y_nn   : y-spline at time step t-1
%             r_nn      : corresponing parameter values at t-1
% 
% OUTPUT      p_nn      : the points at t-1 where the displacement is calculated
%             r_n       : corresponing parameter values at t
%             disp      : the displacment measured in pixels
%             status    : status field indicating success of function call
%                           
% DEPENDENCES   prGetDisplacement uses { SplineLineIntFun
%                                       }
%
% Matthias Machacek 10/09/03

%  _nn variable at timestep t-1
%  _n  variable at timestep t

%calculate the derivative of the last edge spline
dx_sp = fnder(sp_x_nn);
dy_sp = fnder(sp_y_nn);

%determine these locations on the last edge
x_nn=fnval(sp_x_nn,r_nn)';
y_nn=fnval(sp_y_nn,r_nn)';

%determin the normal lines to the last edge:
%  y=ax+b, where a=-dx/dy, b=y-ax
%  the protrusion directions is then -dx,dy
a=(-fnval(dx_sp,r_nn) ./ fnval(dy_sp,r_nn))';
b= y_nn - a .* x_nn;
      
%intersect with actual spline: a x(r)-y(r)+b=0
options = optimset('Display','off');  % Turn off Display
for ii=1:max(size(r_nn))
   r_n(ii) = fsolve(@SplineLineIntFun, r_nn(ii), options, a(ii), b(ii), sp_x_n, sp_y_n);  
end
      
%calculate displacement vectors
x_n=fnval(sp_x_n,r_n)';
y_n=fnval(sp_y_n,r_n)';

p_nn=[x_nn, y_nn];
disp=[x_n-x_nn, y_n-y_nn];

function residual = SplineLineIntFun(x,a,b,sp_x,sp_y)
% SPLINELINEINTFUN calculates the intersection between a spine and a line
% 
% SYNOPSIS    residual = SplineLineIntFun(x,a,b,sp_x,sp_y)
%
% INPUT       x         : parameter of the intersection point
%             a         : parameter defining the line y=ax+b
%             b         : parameter defining the line y=ax+b
%             sp_x      : x-spline structure
%             sp_y      : y-spline structure
% 
% OUTPUT      residual  : the points where the displacement is calculated
%                           
% DEPENDENCES   prGetDisplacement is used by { prGetDisplacement
%                                            }
% Matthias Machacek 10/09/03

 residual = a*fnval(sp_x,x)-fnval(sp_y,x)+b;
