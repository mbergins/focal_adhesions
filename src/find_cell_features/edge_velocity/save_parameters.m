function save_parameters(file_path, parameters)

if isempty(file_path)
    if ispc
        tmp_dir = getenv('TMP'); 
        tmp_dir = [tmp_dir filesep 'parameters.dat'];
    else
        tmp_dir = getenv('HOME');   
        tmp_dir = [tmp_dir filesep '.parameters'];        
    end        
    
    fid = fopen(tmp_dir,'w+');
    if fid == -1
        error('Could not create parameter file');
        return
    end
end

if isfield(parameters,'contr') && ~strcmp(parameters.contr,'empty')
    fprintf(fid, 'CONTR = %f\n',            parameters.contr);
else
    fprintf(fid, 'CONTR = %f\n',            -99);
end 
if isfield(parameters,'file') && ~strcmp(parameters.file,'empty')
    fprintf(fid, 'FILE = %s\n',             parameters.file);
else
    fprintf(fid, 'FILE = %s\n',            'empty');
end   
if isfield(parameters,'results') && ~strcmp(parameters.results,'empty')
    fprintf(fid, 'RESULTS = %s\n',          parameters.results);
else
    fprintf(fid, 'RESULTS = %s\n',          'empty');
end 
if isfield(parameters,'first_img') && ~strcmp(parameters.first_img,'empty')
    fprintf(fid, 'FIRST_IMG = %f\n',        parameters.first_img);
else
    fprintf(fid, 'FIRST_IMG = %f\n',        -99);
end    
if isfield(parameters,'max_img') && ~strcmp(parameters.max_img,'empty')    
    fprintf(fid, 'MAX_IMG = %f\n',          parameters.max_img);
else
    fprintf(fid, 'MAX_IMG = %f\n',          -99);
end  
if isfield(parameters,'t_step') && ~strcmp(parameters.t_step,'empty')
    fprintf(fid, 'T_STEP = %f\n',           parameters.t_step);
else
    fprintf(fid, 'T_STEP = %f\n',           -99);
end    
if isfield(parameters,'movie') && ~strcmp(parameters.movie,'empty')    
    fprintf(fid, 'MOVIE = %f\n',            parameters.movie);
else
    fprintf(fid, 'MOVIE = %f\n',            -99);
end
if isfield(parameters,'bit_depth') && ~strcmp(parameters.bit_depth,'empty')
    fprintf(fid, 'BIT_DEPTH = %f\n',        parameters.bit_depth);
else
    fprintf(fid, 'BIT_DEPTH = %f\n',        -99);
end 
fprintf(fid, 'LAMBDA = %f\n',           0);
fprintf(fid, 'NA = %f\n',               0);
fprintf(fid, 'PIXEL = %f\n',            parameters.pixel);
fprintf(fid, 'TIME_INTERVAL = %f\n',    parameters.time_interval);
if isfield(parameters,'protrusion') && ~strcmp(parameters.protrusion,'empty')
    fprintf(fid, 'PROTRUSION = %f\n',   parameters.protrusion);
else
    fprintf(fid, 'PROTRUSION = %f\n',       -99);
end 
fprintf(fid, 'PROT_SAMPLING = %f\n',    0);
fprintf(fid, 'PROT_DEPTH = %f\n',       0);
fprintf(fid, 'NR_SECT = %f\n',          parameters.nr_sect);
fprintf(fid, 'FILTER_IMAGE = %f\n',     parameters.filter_image);
fprintf(fid, 'IMG_SIGMA = %f\n',        parameters.img_sigma);
fprintf(fid, 'MEDIAN_F = %f\n',         parameters.median_f);
fprintf(fid, 'F_WINDOW = %f\n',         parameters.f_window);
fprintf(fid, 'F_SIGMA = %f\n',          parameters.f_sigma);
fprintf(fid, 'ERODE_DILATE = %f\n',     parameters.erode_dilate);
fprintf(fid, 'TOLERANCE = %f\n',        parameters.tolerance);
fprintf(fid, 'NORMAL = %f\n',           parameters.normal);
fprintf(fid, 'NEAREST = %f\n',          parameters.nearest);
fprintf(fid, 'MECHANICAL = %f\n',       parameters.mechanical);
fprintf(fid, 'TOL = %f\n',              parameters.tol);
fprintf(fid, 'ROBUST_MIN = %f\n',       parameters.robust_min);
fprintf(fid, 'PARENTH_L = %f\n',        parameters.parenth_l);
fprintf(fid, 'PARENTH_R = %f\n',        parameters.parenth_r);
fprintf(fid, 'K_S = %f\n',              parameters.k_s);
fprintf(fid, 'K_W = %f\n',              parameters.k_w);
fprintf(fid, 'CLUSTER = %f\n',          parameters.cluster);
fprintf(fid, 'CLUSTER_METHOD = %s\n',   parameters.cluster_method);
fprintf(fid, 'K_CLUSTER = %f\n',        parameters.k_cluster);
fprintf(fid, 'K_MIN = %f\n',            parameters.k_min);
fprintf(fid, 'K_MAX = %f\n',            parameters.k_min);
fprintf(fid, 'MANUAL_THRESH = %f\n',    parameters.manual_thresh);
fprintf(fid, 'MANUAL_LEVEL = %f\n',     parameters.manual_level);

fclose(fid);