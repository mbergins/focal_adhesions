function prMergePostProc(post_parameters)
% PRMERGEPOSTPROC post processes the results from prMerge    
%
%
% SYNOPSIS          prMergePostProc(post_parameters)
%
% INPUT             varargin  :  varargin       
% 
% OUTPUT               : 
%                           
% DEPENDENCES      prMerge uses {                                
%                                       }
%
%                  prMerge is used by { 
%                                           }
%
% Matthias Machacek 08/01/06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROJECT_DIR          = post_parameters.project_dir;  
MERG_DIR             = post_parameters.merg_dir; 
VAR1_FILE            = post_parameters.var1_file;
VAR2_FILE            = post_parameters.var2_file;
POST_PROC_OPTIONS    = post_parameters.post_proc_options;      
DO_FILTER            = post_parameters.do_filter;
windowSize           = post_parameters.windowSize;
spline_s             = post_parameters.spline_s;
trend                = post_parameters.trend;
corr_lag             = post_parameters.corr_lag;
delay                = post_parameters.delay;


%UNITS = POST_PROC_OPTIONS.units;

FILTER_METHOD = DO_FILTER;
post_parameters.filter_method = FILTER_METHOD;


STATISTICAL_TIME_SERIES_ANALYSIS = 0;
ACTIVITY_MAPS = 0;
CORRELATIONS= 0;
SCATTER_PLOT  = 0;
CURVATURE    = 0;
FFT  = 0;
POWER_SPECTRUM  = 0;
PROT_VS_RET_STAT = 0;
RATIO = 0;

if POST_PROC_OPTIONS == 1
    STATISTICAL_TIME_SERIES_ANALYSIS = 1;
end
if POST_PROC_OPTIONS == 2
    ACTIVITY_MAPS = 1;
end
if POST_PROC_OPTIONS == 3
    CORRELATIONS= 1;
end
if POST_PROC_OPTIONS == 4
    SCATTER_PLOT  = 1;
end
if POST_PROC_OPTIONS == 5
    PROT_VS_RET_STAT = 1;
end
if POST_PROC_OPTIONS == 6
    FFT  = 1;
end
if POST_PROC_OPTIONS == 7
    POWER_SPECTRUM  = 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Important parameters                %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting for the activity maps filtering
post_parameters.img_gauss_w     = 50;
post_parameters.isotropic       = 1;
post_parameters.img_gauss_sig   = 6;
post_parameters.image_stretch   = 15;

DO_ARMA = 0;
%shift the activity, affects CORRELATION
post_parameters.act_time_shift = delay; % positive value shifts forward in time
%shift the scores
post_parameters.sco_time_shift = 0; % positive value shifts forward in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RESULT_DIR  =  MERG_DIR;
post_parameters.merg_dir   =  MERG_DIR;
post_parameters.result_dir =  RESULT_DIR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Get parameters from  merg run  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([MERG_DIR filesep 'merg_parameters.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Get parameters from protrusion run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edge_parameters = read_parameters([post_parameters.project_dir filesep post_parameters.prot_dir filesep 'parameters.dat'], 0);
if edge_parameters.contr == -99
    msgbox('Could not find any edge parameter file','Warning','warn');
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
% determine what units to use
if UNITS == 0
    % use pixel and frame
    PIXEL           = 1;
    TIME_INTERVAL   = 1;
elseif UNITS == 1
    % keep nm and s
    PIXEL           = PIXEL;
    TIME_INTERVAL   = TIME_INTERVAL;
elseif UNITS == 2
    % convert nm to um and sec to min
    PIXEL           =  PIXEL / 1000;
    TIME_INTERVAL   =  TIME_INTERVAL / 60;    
end

if UNITS == 0
    xLabelText = 'frame #';
    yLabelText = 'Velocity (pixel/frame)';
elseif UNITS == 1
    xLabelText = 'Time (s)';    
    yLabelText = 'Velocity (nm/s)';
else  UNITS == 2
    xLabelText = 'Time (min)';    
    yLabelText = 'Velocity (um/min)';
end
end

%x_time_axis = (merg_parameters.first_time: 1: merg_parameters.total_time_steps).* TIME_INTERVAL;


% Create the post processing directory in the merg directory
exist_dir = exist(RESULT_DIR,'dir');
if exist_dir == 0
    %it does not exist, so try to create it
    [s, mess, messid] = mkdir(RESULT_DIR);
    if s==0
        %creation failed, return with error message
        disp('Failed to create the post processing directory');
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load var and var2                                                   %%%%%
%                                                                     %%%%%
% these values are in (pixel/frame) of (nm/s) or (mu/min)             %%%%%
%%%%%% The variables have the structure [segments,time]  %%%%%%%%%%%%%%%%%%

var1_tmp        =  load(VAR1_FILE);
var1_var_name   =  char(fieldnames(var1_tmp));
var1            =  var1_tmp.(var1_var_name);
[pathstr, var1_name, ext, versn] = fileparts(VAR1_FILE);
clear var1_tmp;

var2_tmp        =  load(VAR2_FILE);
var2_var_name   =  char(fieldnames(var2_tmp));
var2            =  var2_tmp.(var2_var_name);
[pathstr, var2_name, ext, versn] = fileparts(VAR2_FILE);
clear var2_tmp;

%get the number of segments
n_segments = size(var1, 1);
%get the number of time steps
n_time     = size(var1, 2);

load([MERG_DIR filesep 'segment_length_av']);
post_parameters.segment_length_av = segment_length_av;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyze statistics of the data                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if STATISTICAL_TIME_SERIES_ANALYSIS
        prStatistical_analysis(edge_parameters, merg_parameters,... 
                        post_parameters, var1, var2, var1_name, var2_name);    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Visualize the data using activity maps                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ACTIVITY_MAPS
    visualize_activity(edge_parameters, merg_parameters, post_parameters,...
                                          var1, var2, var1_name, var2_name);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyze data using FFT algorithm                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FFT
    FFT_activity(edge_parameters, merg_parameters, post_parameters, cell_variables);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Estimate the auto correlations and the                        %
%           cross correlations                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CORRELATIONS
     prCorrelation_analysis(edge_parameters, merg_parameters, post_parameters,...
                                         var1, var2, var1_name, var2_name);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyze the scatter plots                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SCATTER_PLOT
    prScatter_plot(edge_parameters, merg_parameters, post_parameters,...
                                         var1, var2, var1_name, var2_name);       
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Calculates the ratio between activity 1 and 2                 %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RATIO
     ratio_activity(edge_parameters, merg_parameters, post_parameters, cell_variables);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Calculate protrusion vs. retractions statistics               %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PROT_VS_RET_STAT
     prot_vs_ret_activity(edge_parameters, merg_parameters, post_parameters,...
                                         var1, var2, var1_name, var2_name);     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Power spectrum estimators                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if POWER_SPECTRUM
     prPowerSpectrum(edge_parameters, merg_parameters, post_parameters,...
                                         var1, var2, var1_name, var2_name);     
end




