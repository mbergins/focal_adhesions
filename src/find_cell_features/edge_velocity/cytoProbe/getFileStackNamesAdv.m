function [outputFileList]=getFileStackNamesAdv(firstfilename, read_mode)
% getFileStackNames returns a cell array containing all file names (with path) belonging to a stack
%
%
%
%
% SYNOPSIS [outputFileList]=getFileStackNames(firstFileName, read_mode)
%
% INPUT    firstFileName: name of the first greyvalue image to be read
%                    including the full path
%                    the actual filename must consist of 
%                    - alphanumeric body
%                    - numeric number
%                    - extension
%
%           read_mode = 0: read only numeric sequential files
%           read_mode = 1: read all files with numeric part
%           read_mode = 2: read all files
%
% OUTPUT   outputFileList: names of all files belonging to the stack
%                          defined by firstFileName
%
% DEPENDENCES
%
% Aaron Ponti, October 4th, 2002
% Modified by Matthias Machacek 2005

oldDir = [];

% Output
outputFileList = {};

[fpath,fname,fno,fext]=getFilenameBody(firstfilename);

if(~isempty(fname) & ~isempty(fno) & ~isempty(fext) )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % numbered sequence 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(~isempty(fpath))
        % change to stack directory
        oldDir = cd(fpath);
    else
        %check if it is in the matlab search path
        tempName=which(firstfilename);
        if(~isempty(tempName))
            [fpath,fname,fno,fext]=getFilenameBody(tempName);
            oldDir = cd(fpath);
        end
    end

    dirListing = dir;

    % get all relevant filenames
    iEntry = 1;
    fileList = {};
    for( i = 1:length(dirListing))
        if(~dirListing(i).isdir)
            fileList(iEntry) = lower({dirListing(i).name});
            iEntry = iEntry + 1;
        end
    end


    if read_mode == 0
        % get only the numeric sequential file names
        nEntries = 0;
        imIndx = str2num(fno);
        l_fno=length(num2str(fno));
        searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
        outputFileList(1)={strcat(fpath,filesep,searchName)};

        nEntries=1;
        if(~isempty(fileList))
            while( ~isempty(strmatch(lower(searchName),fileList)))
                nEntries = nEntries + 1;
                index(nEntries) = imIndx;
                imIndx = imIndx + 1;
                searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
                outputFileList(nEntries)={strcat(fpath,filesep,searchName)};
            end
        end


        % Removing last file name, which does not exist
        outputFileList=outputFileList(1,1:length(outputFileList)-1);
    elseif read_mode == 1
        % get ALL numeric file names
        nEntries = 0;
        imIndx = str2num(fno);
        l_fno=length(num2str(fno));
        searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
        outputFileList(1)={strcat(fpath,filesep,searchName)};

        nEntries=1;

        % determine last index
        [fpath_last,fname_last,fno_last,fext_last]=getFilenameBody(fileList{end});
        last_index = str2num(fno_last);
        last_index = max(last_index,length(dirListing));

        for i = 1:last_index
            index(nEntries) = imIndx;
            imIndx = imIndx + 1;
            searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
            if(strmatch(lower(searchName),fileList))
                nEntries = nEntries + 1;
                outputFileList(nEntries)={strcat(fpath,filesep,searchName)};
            end
        end
    else
        % get ALL file names
        for ( i = 1:length(dirListing))
            outputFileList = strcat(fpath,filesep, fileList);
        end
    end

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % only one file 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outputFileList{1} = [fpath, filesep, fname,fno,fext];
end

% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end;

