function [outputFileList]=getFileStackNames(firstfilename, exclude_imgs)
% getFileStackNames returns a cell array containing all file names (with path) belonging to a stack
%
% SYNOPSIS [outputFileList]=getFileStackNames(firstFileName)
%
% INPUT    firstFileName: name of the first greyvalue image to be read
%                    including the full path
%                    the actual filename must consist of 
%                    - alphanumeric body
%                    - numeric number
%                    - extension
%          exclude_imgs: optional; vector of image numbers
%                    to exclude
%
% OUTPUT   outputFileList: names of all files belonging to the stack
%                          defined by firstFileName
%
% DEPENDENCES
%
% Aaron Ponti, October 4th, 2002
        
oldDir = [];

% Output
outputFileList = {};

[fpath,fname,fno,fext]=getFilenameBody(firstfilename);

% removed by Matthias to allow for only numeric file names
%if(isempty(fname) | isempty(fno) | isempty(fext) )
%   error('invalid first filename specified');
%end;

% added to allow for individual images stored in 
% numbered directories
indiv_dir_fname = '';
if(isempty(fno))
    [tempPath,tempFname,tempFno,tempFext]=getFilenameBody(fpath);
    if(isempty(tempFno))
        error('invalid first filename specified');
    end
    indiv_dir_fname = strcat(fname,fext);
    fpath = tempPath;
    fname = tempFname;
    fno = tempFno;
    fext = tempFext;
end

if(~isempty(fpath))
	% change to stack directory
   oldDir = cd(fpath);
else
   %check if it is in the matlab search path
   tempName=which(firstfilename);
   if(~isempty(tempName))
      [fpath,fname,fno,fext]=getFilenameBody(tempName);
      
      if(isempty(fno))
        [tempPath,tempFname,tempFno,tempFext]=getFilenameBody(fpath);
        if(isempty(tempFno))
          error('invalid first filename specified');
        end
        indiv_dir_fname = fname + fext;
        fpath = tempPath;
        fname = tempFname;
        fno = tempFno;
        fext = tempFext;
      end
      
      oldDir = cd(fpath);
   end;
end;

dirListing = dir;

% get all relevant filenames
iEntry = 1;
fileList = {};
for(i = 1:length(dirListing))
   if((isempty(indiv_dir_fname) && ~dirListing(i).isdir) || (~isempty(indiv_dir_fname) && dirListing(i).isdir))
      fileList(iEntry) = lower({dirListing(i).name});
      iEntry = iEntry + 1;
   end;
end;

nEntries = 0;
imIndx = str2num(fno);
l_fno=length(num2str(fno));
searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
if(isempty(indiv_dir_fname))
    outputFileList(1)={strcat(fpath,filesep,searchName)};
else
    outputFileList(1)={strcat(fpath,filesep,searchName,filesep,indiv_dir_fname)};
end;
    
nEntries=1;

if(~isempty(fileList))
   while( ~isempty(strmatch(lower(searchName),fileList)))
      nEntries = nEntries + 1;
      index(nEntries) = imIndx;
            
      imIndx = imIndx + 1;
      while (any(imIndx==exclude_imgs))
          imIndx = imIndx + 1;
      end
      
      searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
      if(isempty(indiv_dir_fname))
          outputFileList(nEntries)={strcat(fpath,filesep,searchName)};
      else
          outputFileList(nEntries)={strcat(fpath,filesep,searchName,filesep,indiv_dir_fname)};
      end;
   end;
end;

% Removing last file name, which does not exist 
outputFileList=outputFileList(1,1:length(outputFileList)-1);

% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end;