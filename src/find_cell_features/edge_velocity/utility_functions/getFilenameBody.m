function [path,body,no,ext]=getFilenameBody(fname)
%GETFILENAMEBODY splits a filename into its body, number and 
% extension part
%
% SYNOPSIS [path,body,no,ext]=getFilenameBody(fname)
% 
% INPUT fname : filename; the following filename structure must be 
%                         preserved:
%                         - alphanumeric body
%                         - number before extension
%                         - extension separated
%
% OUTPUT path : string with the path, [] if non-existent
%        body : string with body, [] if non-existent
%        no   : string with the number, [] if non-existent
%        ext  : extension, [] if non-existent
%
% SAMPLE getFileNameBody('test1.tif') returns
%        []
%        'test',
%        '1',
%        '.tif'
%   
%        getFileNameBody('C:\mydir\test1.tif') returns
%        'C:\mydir'
%        'test',
%        '1',
%        '.tif'
%
% SEE ALSO fileparts

% initialize
path = [];
body = [];
no = [];
ext = [];

% search for extension
[path,name,ext,dummy] = fileparts(fname);

% search for letters in remainder
idxs=length(name);
while(idxs > 0 && uint8(name(idxs))>47 && uint8(name(idxs))<58)
   idxs = idxs-1;
end;

% check whether this index points to an actual number
if( length(name) ~= idxs )
   no = name((idxs+1):length(name));
   if(isempty(str2num(no)))
      no = [];
      error('unsupported filename format entered');
   else
      body = name(1:idxs);
   end;
else
   body = name;
end;

   
      
