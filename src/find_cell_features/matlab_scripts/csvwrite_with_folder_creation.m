function csvwrite_with_folder_creation(file_target, values)

[output_folder,~,~] = fileparts(file_target);
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end

csvwrite(file_target,values);

end