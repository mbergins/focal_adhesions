function [p_nn, disp, r_n, nxo_nn, nyo_nn]=prGetDispMechFixL(sp_x_nn, sp_y_nn, sp_x_n, sp_y_n, r_nn, r_n0, contr, varargin)
% PRGETDISPMECH calculates displacement between two two-parametric (x(r),y(r)) splines
% 
%               Given to splines this function calculates the correspondence of 
%               displacement "particles". The correspondence is established
%               based on a mechanical model. The particles at timestep t-1
%               are assumed to be given. The particles at t are connected
%               with each other by a spring. This enforces regular spacing
%               
%
% SYNOPSIS    [p_nn, disp, r_n, nxo_nn, nyo_nn]=prGetDispMechFixL(sp_x_nn, sp_x_nn, sp_x_n, sp_y_n, r_nn, r_n0)
%
% INPUT       sp_x_nn   : x-spline at time step t-1
%             sp_y_nn   : y-spline at time step t-1
%             sp_x_n    : x-spline at time step t
%             sp_y_n    : y-spline at time step t
%             r_nn      : vector with spline parameters where disp. is calc
%             r_n0      : initial protrusion direction parameter values
% 
% OUTPUT      p_nn      : the points at t-1 where the displacement is calculated
%             r_n       : corresponing parameter values at t
%             disp      : the displacment measured in pixels
%             nxo_nn    : the unit oriented (into prot direction) normal 
%             nyo_nn    : the unit oriented (into prot direction) normal
%                           
% DEPENDENCES   prGetDisplacement uses { GlobalMechFun
%                                        GlobalMechFunF
%                                       }
%
%               prGetDisplacement is used by { imEdgeTracker                                     
%                                       }
%
%
% Matthias Machacek 10/24/03

%  _nn variable at timestep t-1
%  _n  variable at timestep t


l=length(varargin);
for i=1:l
    if strcmp(varargin(i),'k_S')
        K_S=varargin{i+1};   
    elseif strcmp(varargin(i),'k_W')
        K_W=varargin{i+1};    
    end
end

%default spring constants
if ~exist('K_S')
    %spacing
    K_S=0.1;
end
if ~exist('K_W')
    %angle
    K_W=1;
end

x_nn=fnval(sp_x_nn,r_nn);
y_nn=fnval(sp_y_nn,r_nn);

%derivative spline 
sp_dx_nn = fnder(sp_x_nn);
sp_dy_nn = fnder(sp_y_nn);
%derivatives at discrite locations
dx_nn=fnval(sp_dx_nn,r_nn);
dy_nn=fnval(sp_dy_nn,r_nn);
%normalize
l=sqrt(dx_nn.^2+dy_nn.^2);
dx_nn=dx_nn./l;
dy_nn=dy_nn./l;
%the normal unit vector (not oriented!!)
nx_nn= dy_nn;
ny_nn=-dx_nn;

%initial values
x0=r_n0;

% check initial values for topological inconsitencies
% [min_val min_index] = min(r_n0);
% [max_val max_index] = max(r_n0);
% if min_index == 1 & max_index == length(r_n0)
%     x0 = sort(r_n0);
%     
% elseif min_index == length(r_n0) & max_index == 1
%     x0 = flipud(sort(r_n0));
%     
% elseif min_index ==  max_index + 1    
%     
%     
% else
%     
%     
% end
    

if contr
    option=optimset('largescale','on','LevenbergMarquardt','on','MaxFunEvals',10000,'TolX',1e-15,'TolFun',1e-15,'display','iter','showstatus','iterplus','gradobj','on','Diagnostics','off');
else
    option=optimset('largescale','on','LevenbergMarquardt','on','MaxFunEvals',10000,'TolX',1e-10,'TolFun',1e-10,'display','off','Diagnostics','off');
end
%[r_n, resnorm, residual_org, exitflag, output] = lsqnonlin(@GlobalMechFun,x0,[],[],option,sp_x_n, sp_y_n, x_nn, y_nn, nx_nn, ny_nn,  K_S, K_W);


if contr
    option2=optimset('largescale','on','MaxFunEvals',10000,'TolX',1e-15,'TolFun',1e-15,'display','iter','showstatus','iterplus','gradobj','on','Diagnostics','off');
else
    option2=optimset('largescale','on','MaxFunEvals',10000,'TolX',1e-15,'TolFun',1e-15,'display','off','Diagnostics','off');
end
[r_n, fval, exitflag, output] = fsolve(@GlobalMechFun,x0,option2,sp_x_n, sp_y_n, x_nn, y_nn, nx_nn, ny_nn,  K_S, K_W);

% here also the non-oriented normal x_nn is oriented: nxo_nn
[F, S1_t, S2_t, M, W_t, nxo_nn, nyo_nn] = GlobalMechFunF(r_n, sp_x_n, sp_y_n, x_nn, y_nn, nx_nn, ny_nn, K_S, K_W);
S=S1_t-S2_t;
W=W_t;

residual=F;


if contr
	figure;
	subplot(3,1,1);
    plot(residual);
	title('Residual: force on protrusion particle');
	text(2,0.8*max(residual),['Weights : Spacing, K_S=  ',num2str(K_S)],'Color','r');
	text(2,0.4*max(residual),['Angle, K_S=  ',num2str(K_W)],'Color','r');
	subplot(3,1,2);
    plot(S);
	title('Distance force');
	subplot(3,1,3);
    plot(W);
	title('Angular force');
end

%calculate displacement vectors
x_n=fnval(sp_x_n,r_n)';
y_n=fnval(sp_y_n,r_n)';

p_nn=[x_nn', y_nn'];
disp=[x_n-x_nn', y_n-y_nn'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=GlobalMechFun(r_n, s_x_n, s_y_n, x_nn, y_nn, nx_nn, ny_nn, K_S, K_W)

%  _nn variable at timestep t-1
%  _n  variable at timestep t
knots_nr=length(r_n);

% get the length of the spline covered by the knots
total_spring_length = max(r_n)-1;
spring_length = total_spring_length/length(r_n)-1;

% positions at the new spline
x_n=fnval(s_x_n,r_n);
y_n=fnval(s_y_n,r_n);

% derivative new spline 
s_dx_n = fnder(s_x_n);
s_dy_n = fnder(s_y_n);

% derivatives at the new spline
dx_n=fnval(s_dx_n,r_n);
dy_n=fnval(s_dy_n,r_n);

% spacing vector to the right particle
for i=1:knots_nr-1
   l1_x(i)= x_n(i+1) - x_n(i);
   l1_y(i)= y_n(i+1) - y_n(i);
end

% spacing vector length
l=sqrt(l1_x.^2+l1_y.^2);  

%if there is zero element, set the value to 0.0001
%this might happen when the spline parameter iterates out
%of its range
[i_l, j_l, val_l] = find(l==0);
l(j_l) = 0.00001; 

%unit spacing vector to the right particle
l1_x_u=l1_x./l;
l1_y_u=l1_y./l;

%the deflection from its neutral position: delta_x
for i=1:knots_nr-1
   delta_x(i)= spring_length - l(i);
end


if K_S ~= 0
    %force from the spring elements
    for i=1:knots_nr-1
%         if r_n(i) > r_n(i+1)
%             S_x(i)= 10000 * K_S  * abs(delta_x(i)^2)  * l1_x_u(i);
%             S_y(i)= 10000 * K_S  * abs(delta_x(i)^2)  * l1_y_u(i);                
%         else
%             S_x(i)= K_S  * delta_x(i)  * l1_x_u(i);
%             S_y(i)= K_S  * delta_x(i)  * l1_y_u(i);
%         end    
        S_x(i)= K_S  * delta_x(i)  * l1_x_u(i);
        S_y(i)= K_S  * delta_x(i)  * l1_y_u(i);
    end
    
    % get the normalized tangent at the new edge 
    ld=sqrt(dx_n.^2+dy_n.^2);

    %if there is zero element, set the value to 0.0001
    %this might happen when the spline parameter iterates out
    %of its range
    [i_l, j_l, val_l] = find(ld==0);
    ld(j_l) = 0.00001; 
    
    dx_n=dx_n./ld;
    dy_n=dy_n./ld;
      
    for i=2:knots_nr-1
        % force from the left spring
        S1_t(i) = S_x(i-1) * dx_n(i) + S_y(i-1) * dy_n(i);
        % force from the right spring
        S2_t(i) = S_x(i)   * dx_n(i) + S_y(i)   * dy_n(i); 
    end
    
    % boundery conditions
    S1_t(1) = S1_t(2);
    S2_t(1) = S_x(1) * dx_n(1) + S_y(1) * dy_n(1);
    S1_t(knots_nr) = S_x(knots_nr-1) * dx_n(knots_nr) + S_y(knots_nr-1) * dy_n(knots_nr);
    S2_t(knots_nr) = S2_t(knots_nr-1);
    
else 
    S1_t = 0;
    S2_t = 0;
end

if K_W ~= 0
    %distance from the particle at the last time step
    dist = sqrt((x_nn - x_n).^2 + (y_nn - y_n).^2);
    % 
    %%%%%%%%%  force from the angular deviation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit protrusion direction
    d_x = x_n-x_nn;
    d_y = y_n-y_nn;
    d=sqrt(d_x.^2+d_y.^2);
    d_x = d_x ./ d;
    d_y = d_y ./ d;
    %first determine the direction of the normal vector 
    %so that it is pointing in the protrusion direction
    sg=sign(d_x.*nx_nn+d_y.*ny_nn);
    for i=1:knots_nr
        if sg(i)<0
            nx_nn(i)=-nx_nn(i);
            ny_nn(i)=-ny_nn(i);
        end
    end 
    %sign of the momentum force : sign(a x b)
    sg=sign(nx_nn.*d_y - ny_nn.*d_x);
    %M=sg.*0.5.*(ones(1,length(sg))-(nx_nn.*d_x + ny_nn.*d_y));
    angle_value = nx_nn.*d_x + ny_nn.*d_y;
    % check for values bigger than 1
    [i_angle_value, j_angle_value, val_angle_value] = find(abs(angle_value) > 1);
    angle_value(j_angle_value) = 1; 
    M = sg./pi.* acos(angle_value);
    %direction of the angular force
    w_x =  d_y;
    w_y = -d_x;
    %magnitude of the angular force, here (tx_nn, ty_nn) is the tangent 
    W_x =  K_W .* M .* w_x;
    W_y =  K_W .* M .* w_y;
    
    % project the angluar force onto the tangent
    W_t = W_x .* dx_n + W_y .* dy_n;
else
    M   = 0;
    W_t = 0;
end


F   =   S1_t - S2_t + W_t;

if sum(isnan(F))
   questdlg('NaN in mechanic model edge tracker. Continue?','Warning');
end
if sum(isinf(F))
    questdlg('Inf in mechanic model edge tracker. Continue?','Warning');
end
    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, S1_t, S2_t, M, W_t, nx_nn, ny_nn] = GlobalMechFunF(r_n, s_x_n, s_y_n, x_nn, y_nn, nx_nn, ny_nn, K_S, K_W)

%  _nn variable at timestep t-1
%  _n  variable at timestep t
knots_nr=length(r_n);

% get the length of the spline covered by the knots
total_spring_length = max(r_n)-1;
spring_length = total_spring_length/length(r_n)-1;

x_n=fnval(s_x_n,r_n);
y_n=fnval(s_y_n,r_n);


% derivative new spline 
s_dx_n = fnder(s_x_n);
s_dy_n = fnder(s_y_n);

% derivatives at the new spline
dx_n=fnval(s_dx_n,r_n);
dy_n=fnval(s_dy_n,r_n);

if K_S ~= 0
    %spacing vector to the right particle
    for i=1:knots_nr-1
        l1_x(i)= x_n(i+1) - x_n(i);
        l1_y(i)= y_n(i+1) - y_n(i);
    end
    %spacing vector length
    l=sqrt(l1_x.^2+l1_y.^2);   
    
    %if there is zero element, set the value to 0.0001
    %this might happen when the spline parameter iterates out
    %of its range
    [i_l, j_l, val_l] = find(l==0);
    l(j_l) = 0.00001; 
    
    %unit spacing vector to the right particle
    l1_x_u=l1_x./l;
    l1_y_u=l1_y./l;
    
    %the deflection from its neutral position: delta_x
    for i=1:knots_nr-1
        delta_x(i)= spring_length - l(i);
    end
    
    %force from the spring elements
    for i=1:knots_nr-1
        S_x(i)= K_S * delta_x(i) * l1_x_u(i);
        S_y(i)= K_S * delta_x(i) * l1_y_u(i);
    end
    %S_x(knots_nr) = S_x(knots_nr-1);
    %S_y(knots_nr) = S_y(knots_nr-1);
    
    % get the normalized tangent at the new edge 
    ld=sqrt(dx_n.^2+dy_n.^2);
    dx_n=dx_n./ld;
    dy_n=dy_n./ld;
   
    
    for i=2:knots_nr-1
        % force from the left spring
        S1_t(i) = S_x(i-1) * dx_n(i) + S_y(i-1) * dy_n(i);
        % force from the right spring
        S2_t(i) = S_x(i) *   dx_n(i) + S_y(i) *   dy_n(i); 
    end
    S1_t(1) = S1_t(2);
    S2_t(1) = S_x(1)   * dx_n(1) + S_y(1)   * dy_n(1);
    S1_t(knots_nr) = S_x(knots_nr-1) * dx_n(knots_nr) + S_y(knots_nr-1) * dy_n(knots_nr);
    S2_t(knots_nr) = S2_t(knots_nr-1);
else
    S1_t = 0;
    S2_t = 0;    
end


if K_W ~= 0
   %distance from the particle at the last time step
    dist = sqrt((x_nn - x_n).^2 + (y_nn - y_n).^2);    
    %%%%%%%%%%%%%% force from the angular deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit protrusion direction
    d_x = x_n-x_nn;
    d_y = y_n-y_nn;
    d=sqrt(d_x.^2+d_y.^2);
    d_x = d_x ./ d;
    d_y = d_y ./ d;
    % first determine the direction of the normal vector 
    % so that it is pointing in the protrusion direction
    sg=sign(d_x.*nx_nn+d_y.*ny_nn);
    for i=1:knots_nr
        if sg(i)<0
            nx_nn(i)=-nx_nn(i);
            ny_nn(i)=-ny_nn(i);
        end
    end 
    % sign of the momentum force : sign(a x b)
    sg=sign(nx_nn.*d_y - ny_nn.*d_x);
    % M=sg.*0.5.*(ones(1,length(sg))-(nx_nn.*d_x + ny_nn.*d_y));
    angle_value = nx_nn.*d_x + ny_nn.*d_y;
    [i_angle_value, j_angle_value, val_angle_value] = find(abs(angle_value) > 1);
    angle_value(j_angle_value) = 1; 
    M= sg./pi.*acos(angle_value);
    % direction of the angular force
    w_x =  d_y;
    w_y = -d_x;
    % magnitude of the angular force, here (tx_nn, ty_nn) is the tangent 
    W_x = K_W .* M .* w_x;
    W_y = K_W .* M .* w_y;
    
    % project the angluar force onto the tangent
    W_t = W_x .* dx_n + W_y .* dy_n;
else
    M   = 0;
    W_t = 0;
end


F   =   S1_t - S2_t + W_t;
