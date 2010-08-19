function save_pralpha_parameters(file_path, parameters)

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
else
    fid = fopen(file_path,'w+');
    if fid == -1
        error('Could not create parameter file');
        return
    end       
end

fprintf(fid, 'PROJECT_DIR           = %s\n',    parameters.project_dir);
fprintf(fid, 'IMG_NAME              = %s\n',    parameters.img_name);
fprintf(fid, 'PROT_DIR              = %s\n',    parameters.prot_dir);
fprintf(fid, 'FLOW_DIR              = %s\n',    parameters.flow_dir);
fprintf(fid, 'SCORES_DIR            = %s\n',    parameters.scores_dir);
fprintf(fid, 'ACTIVITY_DIR_1        = %s\n',    parameters.activity_dir_1);
fprintf(fid, 'ACTIVITY_DIR_2        = %s\n',    parameters.activity_dir_2);
fprintf(fid, 'DO_FLOW               = %i\n',    parameters.do_flow);
fprintf(fid, 'DO_SCORES             = %i\n',    parameters.do_scores);
fprintf(fid, 'DO_PROT               = %i\n',    parameters.do_prot);
fprintf(fid, 'DO_ACTIVITY_1         = %i\n',    parameters.do_activity_1);
fprintf(fid, 'DO_ACTIVITY_2         = %i\n',    parameters.do_activity_2);
fprintf(fid, 'FIRST_TIME            = %i\n',    parameters.first_time);
fprintf(fid, 'TOTAL_TIME_STEPS      = %i\n',    parameters.total_time_steps);
fprintf(fid, 'START_SEG             = %i\n',    parameters.start_seg);
fprintf(fid, 'END_SEG               = %i\n',    parameters.end_seg);
fprintf(fid, 'SEG_SHIFT             = %f\n',    parameters.seg_shift);
fprintf(fid, 'INTERPOLATE           = %i\n',    parameters.interpolation);
fprintf(fid, 'SEG_NR                = %i\n',    parameters.seg_nr);
fprintf(fid, 'SEG_LENGTH            = %f\n',    parameters.seg_length);
fprintf(fid, 'SEG_DEPTH             = %f\n',    parameters.seg_depth);
fprintf(fid, 'D0                    = %f\n',     parameters.d0);
fclose(fid);
