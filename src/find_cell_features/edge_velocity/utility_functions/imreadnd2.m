function I=imreadnd2(filename,xmin,xmax)
% IMREADND2 reads an image file and normalizes all pixel intensities in the range [xmin..xmax] -> [0..1]
%
% SYNOPSIS   I=imreadnd2(filename,xmin,xmax)
%
% INPUT      filename : name of the file with entire path
%            xmin     : minimum gray value used to stretch the intensities (e.g. 0)
%            xmax     : maximum gray value used to stretch the intensities (e.g. 255)
%
% OUTPUT     I        : image with intensities stretched to the [0..1] range
%
% Aaron Ponti, 2001

if nargin<3
   error('This functions requires three parameter');   
end

% 'filename' is passed 'as is' to imread, which will check it
X=imread(filename);

% Conversion of I to type double (necessary for normalization)
X=double(X);

% Normalization
I=(X-xmin)/(xmax-xmin);

% Check for normalization
ch=max(I(:));
if ch>1
   error('ERROR: Normalization failed! Select a greater bit depth.');
end

