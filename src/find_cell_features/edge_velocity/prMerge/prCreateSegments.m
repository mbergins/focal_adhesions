function [xb, yb] = prCreateSegments(sp_x, sp_y, x_n, y_n, m_img, n_img, SEG_NR, SEG_DEPTH, i_seg)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Create sectors along the edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



max_spline_par = sp_x.knots(end);
% s_p = (max_spline_par-1)/SEG_NR;
% m_pos=1:s_p:max_spline_par;

%Customized by Lin Ji on Sept. 21 to make the cutoff segments' width at
% the beginning and end smaller.
s_p = (max_spline_par-1)/SEG_NR;
cutOffSegWidth = min(10,max(1,s_p/2));
m_pos = linspace(cutOffSegWidth+1,max_spline_par-cutOffSegWidth,SEG_NR+1);

% extract image segment
spl=m_pos(i_seg):m_pos(i_seg+1);
xb=fnval(sp_x,spl);
yb=fnval(sp_y,spl);
le=length(xb);
for i2=1:le
    xb(le+i2) = xb(le+1-i2) + SEG_DEPTH*(-x_n(i_seg));
    yb(le+i2) = yb(le+1-i2) + SEG_DEPTH*(-y_n(i_seg));
    if xb(le+i2) < 1
        xb(le+i2) = 1;
    elseif xb(le+i2) > m_img
        xb(le+i2) = m_img;
    end
    if yb(le+i2) < 1
        yb(le+i2) = 1;
    elseif yb(le+i2) > n_img
        yb(le+i2) = n_img;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%