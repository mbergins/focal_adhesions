function dist = distanceFromLineSegment(x,y,p0,p1)
% dist = distanceFromLineSegment(x,y,p0,p1)
% dist(i,j) = distance of [x(i,j),y(i,j)] to line segment [p0,p1]
% p0,p1 are 2x1 vectors, and are of the form (x,y)
% Example:
%  [x y] = meshgrid(1:360,1:240);
%  p0    = [10,20]; %(x,y)
%  p1    = [300,120]; %(x,y)
%  dist  = distanceFromLineSegment(x,y,p0,p1);
%  %draw a (anti-aliased) line corresponding to the segment
%  %exp(-r^2/(2*sigma^2)),
%  line = exp(-dist.^2 / 16);
%  subplot(121); imagesc(dist); axis equal;
%  subplot(122); imshow(line);
%

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

% orthogonal distance from the line segment (as per MathWorld)
% http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    % for easier typing and reading
    x2 = p1(1); y2 = p1(2);
    x1 = p0(1); y1 = p0(2);

    % equation of the plane
    % ax + by + c = 0
    % therefore, c = -by - ax
    b = -(x2-x1); a =  (y2-y1); c =  -b*y1 - a*x1;

    numer = abs(a*x + b*y + c);
    denom = sqrt((x2-x1).^2 + (y2-y1).^2);
    orthoDist = numer ./ denom;

% Distance from point p0 and p1
    distFromP0 = sqrt((x-x1).^2 + (y-y1).^2);
    distFromP1 = sqrt((x-x2).^2 + (y-y2).^2);

% figure out where to apply each distance appropriately
    % perpLine1 a1x + b1y + c1 = 0
    b1 = a; a1 = -b; c1 = -b1*y1 - a1*x1;
    % which direction is determined by testing the other point
    correctSign = -sign(a1*x2 + b1*y2 + c1);
    mask1 = sign(a1*x + b1*y + c1) == correctSign;
    %similarly, perpLine2 a2x + b2y + c2 = 0
    %(same for a1=a2,b1=b2, but repeated for clarity)
    b2 = a; a2 = -b; c2 = -b2*y2 - a2*x2;
    % which direction is determined by testing the other point
    correctSign = -sign(a2*x1 + b2*y1 + c2);
    mask2 = sign(a2*x + b2*y + c2) == correctSign;

% report the distance based on the coordinates
    dist = orthoDist;
    dist(mask1) = distFromP0(mask1);
    dist(mask2) = distFromP1(mask2);
