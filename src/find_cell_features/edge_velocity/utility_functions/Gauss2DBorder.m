function [out_corrected,M]=Gauss2DBorder(x,sigma);
% Gauss2DBorder	apply a 2 dimensional gauss filter on enlarged image
%
%    [out_corrected,M]=Gauss2DBorder(x,sigma);
%
%    INPUT: x                   image
%           sigma               of gauss filter
%
%    OUTPUT: out_corrected      filtered image
%            M                  gaussian mask
%

% bug fix: AP - 10.07.02
% Added correction for image boarder effects
% Matthias Machacek 04/05/04

R = ceil(3*sigma);   % cutoff radius of the gaussian kernel 

for i = -R:R,
   for j = -R:R,
      M(i+R+1,j+R+1) = exp(-(i*i+j*j)/2/sigma/sigma);
   end
end
M = M/sum(M(:));   % normalize the gaussian mask so that the sum is
                   % equal to 1

 
% initialize new image
[img_h img_w] = size(x);
x_extendet = zeros(img_h + 2*R, img_w + 2*R);

% copy image 
x_extendet(R+1 : end-R,  R+1 : end-R) = x;

% left boarder
x_extendet(R+1 : end-R,1:R)             = fliplr(x(:,1:R));
% right boarder
x_extendet(R+1 : end-R,end-R:end)       = fliplr(x(:,end-R:end));
% upper boarder
x_extendet(1 : R, R+1 : end-R)          = flipud(x(1 : R,:));
% lower boarder
x_extendet(end-R : end, R+1 : end-R)    = flipud(x(end-R : end,:));

% fill upper left corner
x_extendet(1:R ,1:R)                    = fliplr(flipud(x(1 : R,1:R)));
% fill upper right corner
x_extendet(1:R ,end-R:end)              = fliplr(flipud(x(1 : R,end-R:end)));
% fill lower left corner
x_extendet(end-R:end ,1:R)              = fliplr(flipud(x(end-R:end ,1:R)));
% fill lower right corner
x_extendet(end-R:end ,end-R:end)        = fliplr(flipud(x(end-R:end ,end-R:end)));


% Convolute matrices
out=filter2(M,x_extendet);
out_corrected=out(R+1 : end-R,  R+1 : end-R);
