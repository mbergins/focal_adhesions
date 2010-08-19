function [A_av,A_std] =  block_average(A,dt,av)

% A(seg,time)
% dt: time step
% av: average window

T = size(A,2)-av;
j=1;

for i1=1:dt:T
    i2 = i1+av-1;
    
    A_av(:,j)  = sum(A(:,i1:i2),2) / av;
    A_std(:,j) = std(A(:,i1:i2),0,2);
    
    j=j+1;
end

protrusion_normal = A_av;
save('A_av','protrusion_normal');
protrusion_normal = A_std;
save('A_std','protrusion_normal');

% save
