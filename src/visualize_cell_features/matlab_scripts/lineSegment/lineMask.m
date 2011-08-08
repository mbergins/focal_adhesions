function mask = lineMask(imgSize, p0, p1, sigma)
% function mask = lineMask(imgSize, p0, p1 [, sigma])
% generates a weighted line mask(i,j) = weight,
% optional argument sigma (default=4) controls the width of the line
% imgSize = [height, width]
% p0 = [row,col]
% p1 = [row,col] 
% points do not have to be within the image
% example:
%  img  = rgb2gray(im2double(imread('peppers.png')));
%  mask = lineMask(size(img), [50,50], [1005,4000]);
%  %put a gray (intensity=.75) line in the image at the place specified
%  out  = img.*(1-mask) + .75*mask;
%  figure(1); imshow(mask);
%  figure(2); imshow(out);
% REQUIRES distanceFromLineSegment.m

%
% by alvin.a.raj@gmail.com 11/23/2009

% Copyright (c) 2009 Alvin Raj
% 
%  Permission is hereby granted, free of charge, to any person
%  obtaining a copy of this software and associated documentation
%  files (the "Software"), to deal in the Software without
%  restriction, including without limitation the rights to use,
%  copy, modify, merge, publish, distribute, sublicense, and/or sell
%  copies of the Software, and to permit persons to whom the
%  Software is furnished to do so, subject to the following
%  conditions:
% 
%  The above copyright notice and this permission notice shall be
%  included in all copies or substantial portions of the Software.
% 
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
%  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%  OTHER DEALINGS IN THE SOFTWARE.


if (~exist('sigma','var'))
    sigma = 4;
end

[x y] = meshgrid(1:imgSize(2),1:imgSize(1));
% convert (row,col) -> (x,y) 
p0 = [p0(2),p0(1)]; 
p1 = [p1(2),p1(1)]; 
dist = distanceFromLineSegment(x,y,p0,p1);
mask = exp(-dist.^2 / (2*sigma^2));