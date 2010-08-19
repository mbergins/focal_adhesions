function parameters = read_parameters(file_path, write_tmp)

WRITE_TMP = write_tmp;

if WRITE_TMP
    if ispc
        tmp_dir = getenv('TMP'); 
        tmp_dir = [tmp_dir filesep 'parameters.dat'];        
    else
        tmp_dir = getenv('HOME');   
        tmp_dir = [tmp_dir filesep '.parameters'];        
    end

    
    fid = fopen(tmp_dir,'r');
    if fid == -1
        parameters.contr = -99;
        return
    end
else
    fid = fopen(file_path,'r');
    if fid == -1
        parameters.contr = -99;
        return
    end
end


s = fscanf(fid, '%s', [1 2]);
parameters.contr                =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.file                 =   fscanf(fid, '%s', [1 1]);
s = fscanf(fid, '%s', [1 2]);
parameters.results              =   fscanf(fid, '%s', [1 1]);
s = fscanf(fid, '%s', [1 2]);
parameters.first_img            =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.max_img              =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.t_step               =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.movie                =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.bit_depth            =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.lambda               =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.na                   =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.pixel                =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.time_interval        =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.protrusion           =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.prot_sampling        =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.prot_depth           =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.nr_sect              =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.filter_image         =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.img_sigma            =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.median_f             =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.f_window             =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.f_sigma              =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.erode_dilate         =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.tolerance            =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.normal               =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.nearest              =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.mechanical           =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.tol                  =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.robust_min           =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.parenth_l            =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.parenth_r            =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.k_s                  =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.k_w                  =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.cluster              =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.cluster_method       =   fscanf(fid, '%s', [1 1]);
s = fscanf(fid, '%s', [1 2]);
parameters.k_cluster            =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.k_max                =   str2num(fscanf(fid, '%s', [1 1]));   
s = fscanf(fid, '%s', [1 2]);
parameters.k_min                =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.manual_thresh        =   str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.manual_level         =   str2num(fscanf(fid, '%s', [1 1]));
fclose(fid);