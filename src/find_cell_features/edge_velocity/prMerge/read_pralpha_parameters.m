function parameters = read_pralpha_parameters(file_path)

fid = fopen(file_path,'r');
if fid == -1
    parameters.first_time = -1;
    return
end

s = fscanf(fid, '%s', [1 2]);
parameters.prot_dir = fscanf(fid, '%s', [1 1]);
s = fscanf(fid, '%s', [1 2]);
parameters.img_name = fscanf(fid, '%s', [1 1]);
s = fscanf(fid, '%s', [1 2]);
parameters.prot_dir = fscanf(fid, '%s', [1 1]);
s = fscanf(fid, '%s', [1 2]);
parameters.fsm_dir = fscanf(fid, '%s', [1 1]);
s = fscanf(fid, '%s', [1 2]);
parameters.first_time = str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.total_time_steps = str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.start_seg = str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.end_seg = str2num(fscanf(fid, '%s', [1 1]));
s = fscanf(fid, '%s', [1 2]);
parameters.seg_shift = str2num(fscanf(fid, '%s', [1 1]));



% save file info
% fprintf(fid, 'PROJECT_DIR           = s\n',    parameters.project_dir);
% fprintf(fid, 'IMG_NAME              = s\n',    parameters.img_name);
% fprintf(fid, 'PROT_DIR              = s\n',    parameters.prot_dir);
% fprintf(fid, 'FSM_DIR               = s\n',    parameters.fsm_dir);
% fprintf(fid, 'FIRST_TIME            = f\n',    parameters.first_time);
% fprintf(fid, 'TOTAL_TIME_STEPS      = f\n',    parameters.total_time_steps);
% fprintf(fid, 'START_SEG             = f\n',    parameters.start_seg);
% fprintf(fid, 'END_SEG               = f\n',    parameters.end_seg);
% fprintf(fid, 'SEG_SHIFT             = f\n',    parameters.seg_shift);