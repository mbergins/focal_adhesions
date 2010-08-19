function [p_nn, disp, r_n, status]=prGetDispNearest(sp_x_nn, sp_y_nn, sp_x_n, sp_y_n, r_nn, varargin)
% PRGETDISPLACEMENT calculates displacement between two two-parametric (x(r),y(r)) splines
%
%
%             If you want to be really sure that you find the shortest
%             displacement set ROBUST_MIN=1, and TOL=0. In this case the
%             entire spline will be searched for the shortest displacement.
%             This might be slow. 
%             You can limit the search space by setting TOL ~= 0. In this
%             case the global minima is found in a range [r_nn(ii) - TOL, 
%             r_nn(ii) + TOL], where r_nn(ii) is the spline parameter in
%             the last time step. Be aware that the more unrelated to
%             spline parameter between two timesteps is the bigger TOL has
%             to be. 
%             In case of ROBUST_MIN=1 the accuracy of the solution is given
%             by SOLUTION_SPACING.
%             For ROBUST_MIN the minimal distance is found by the matlab
%             function fminbnd. The solution can be a local solution!!!
% 
% SYNOPSIS    [p_nn, disp, r_n, status]=prGetDispNearest(sp_x_nn, sp_x_nn, sp_x_n, sp_y_n, r_nn)
%
% INPUT       sp_x_n    : x-spline at time step t+1
%             sp_y_n    : y-spline at time step t+1
%             sp_x_nn   : x-spline at time step t
%             sp_y_nn   : y-spline at time step t
%             r_nn      : corresponing parameter values at t
% 
% OUTPUT      p_nn      : the points at t where the displacement is calculated
%             r_n       : corresponing parameter values at t+1
%             disp      : the displacment measured in pixels
%             status    : status field indicating success of function call
%                           
% DEPENDENCES   prGetDisplacement uses { SplineLineDistFun
%                                       }
%
% Matthias Machacek 10/11/03

l=length(varargin);
for i=1:2:l
    in_found=0;
    if strcmp(varargin(i),'tol')
        TOL=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'robust_min')
        ROBUST_MIN=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'solution_spacing')
        SOLUTION_SPACING=varargin{i+1};
        in_found=1;         
    end
    
    if in_found == 0
        error_string = char(varargin(i));
        error(['Unknown input:   ' , error_string]);
    end       
end

if ~exist('TOL','var')
    TOL=12;
end
if ~exist('ROBUST_MIN','var')
    ROBUST_MIN=0;
end
if ~exist('SOLUTION_SPACING','var')
    SOLUTION_SPACING = 0.01;
end


if ROBUST_MIN
    if TOL == 0
        r_n_lower = sp_x_n.knots(1);
        r_n_upper = sp_x_n.knots(end);      
        r_n_cand = r_n_lower: SOLUTION_SPACING : r_n_upper;
    end
    
    for ii=1:size(r_nn,2)
        %determine these locations on the last edge
        x_nn(ii) = fnval(sp_x_nn,r_nn(ii));
        y_nn(ii) = fnval(sp_y_nn,r_nn(ii));
        
        if TOL ~= 0
            %the lower boundery of the parameter search space
            r_n_lower = r_nn(ii) - TOL; 
            if r_n_lower < sp_x_n.knots(1)
                r_n_lower = sp_x_n.knots(1);
            end
            %the upper boundery of the parameter search space
            r_n_upper = r_nn(ii) + TOL; 
            if r_n_upper > sp_x_n.knots(end)
                r_n_upper = sp_x_n.knots(end);
            end
            r_n_cand = r_n_lower: SOLUTION_SPACING : r_n_upper;
        end
        
        
%         %for the first point consider the whole edge
%         if ii == 1
            %create possible candidates
            
%         else
%             %now create possible candidates based on last found candidate
%             %the lower boundery of the parameter search space
%             r_n_lower = r_n(ii-1) - TOL; 
%             if r_n_lower < sp_x_nn.knots(1)
%                 r_n_lower = sp_x_nn.knots(1);
%             end
%             %the upper boundery of the parameter search space
%             r_n_upper = r_n(ii-1) + TOL; 
%             if r_n_upper > sp_x_nn.knots(end)
%                 r_n_upper = sp_x_nn.knots(end);
%             end
%             
%             r_n_cand = r_n(ii-1)-r_n_lower: SOLUTION_SPACING : r_n(ii-1)+r_n_upper;
%             
%         end
 
        if ii   ==  size(r_nn,2)-3
            a=1;
        end
            
        %calculate distances
        dist_cand = sqrt((fnval(sp_x_n, r_n_cand)- x_nn(ii)).^2+ (fnval(sp_y_n, r_n_cand)-y_nn(ii)).^2);
 
        %get the minimal distance
        [dist min_index] = min(dist_cand);
        
        if length(min_index) == 1
            %get the corresponding parameter
            r_n(ii) = r_n_cand(min_index);    
        else
            a = 1;
        end
        
        %test if the minimal distance is in the parameter range
        if r_n(ii) < 1 | r_n(ii) > sp_x_n.knots(end)
            return
        end
    end
else
    options = optimset('Display','off');  % 'notify' Turn off Display
    options = [];
    for ii=1:size(r_nn,2)
        %determine these locations on the last edge
        x_nn(ii) = fnval(sp_x_nn,r_nn(ii));
        y_nn(ii) = fnval(sp_y_nn,r_nn(ii));
        
        
        %test the boudery limits
        lower_limit = r_nn(ii)-TOL;
        if lower_limit < sp_x_nn.knots(1)
            lower_limit = sp_x_nn.knots(1);
        end
        upper_limit = r_nn(ii)+TOL;
        if upper_limit > sp_x_nn.knots(end)
            upper_limit = sp_x_nn.knots(end);
        end
        
        [r_n(ii), fval, status] = fminbnd(@SplineLineDistFun, lower_limit, upper_limit, options, sp_x_n, sp_y_n, x_nn(ii), y_nn(ii));
        %r_n(ii) = fsolve(@SplineLineDistFun, r_nn(ii), options, sp_x_n, sp_y_n, x_nn, y_nn);
        
        
        %check if a solution was found
        %    if fval < 0
        %         r_n(ii)  = -99;  
        %    end
        %    if fval == 0
        %         r_n(ii)  = -99;  
        %    end  
    end
end


%calculate displacement vectors
x_n=fnval(sp_x_n,r_n)';
y_n=fnval(sp_y_n,r_n)';

p_nn=[x_nn', y_nn'];
disp=[x_n-x_nn', y_n-y_nn'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist=SplineLineDistFun(r_n, sp_x_n, sp_y_n, x_nn, y_nn)

dist=sqrt((fnval(sp_x_n,r_n)-x_nn)^2+(fnval(sp_y_n,r_n)-y_nn)^2);



