function max_image_num = find_max_image_num(folder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('folder',@(x)exist(x,'dir')==7);

i_p.parse(folder);

folders_to_exclude = {'.','..'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder_dir = dir(folder);
max_image_num = -Inf;
for i = 1:size(folder_dir)
    if (strmatch(folder_dir(i).name,folders_to_exclude))
        continue;
    end
    
    if (str2num(folder_dir(i).name) > max_image_num), max_image_num = str2num(folder_dir(i).name); end
end