function varargout = prMergePanel(varargin)
% PRMERGEPANEL M-file for prMergePanel.fig
%      PRMERGEPANEL, by itself, creates a new PRMERGEPANEL or raises the existing
%      singleton*.
%
%      H = PRMERGEPANEL returns the handle to a new PRMERGEPANEL or the handle to
%      the existing singleton*.
%
%      PRMERGEPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRMERGEPANEL.M with the given input arguments.
%
%      PRMERGEPANEL('Property','Value',...) creates a new PRMERGEPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prMergePanel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prMergePanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prMergePanel

% Last Modified by GUIDE v2.5 10-Apr-2007 10:26:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prMergePanel_OpeningFcn, ...
                   'gui_OutputFcn',  @prMergePanel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before prMergePanel is made visible.
function prMergePanel_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.pr_alpha_option(1),'Check','on');

%fill the text boxes 
handles.parameter_set = 0;
set(handles.total_time_steps,'String',10);
set(handles.first_time,'String',1);
handles.debug = 0;

% fill the segment definition list
vars{1} = 'Number of segments';
vars{2} = 'Segment length [pixel]';
set(handles.seg_def,'String',vars);
set(handles.seg_val,'String',30);
set(handles.seg_depth,'String',10);

% Influence diameter for the flow field
% averaging from the FSM data
set(handles.d0,'String',10);
set(handles.d0,'Enable','off');

%%%%%%%%%% Postprocessing elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.windowSize,'String',20);
set(handles.spline_s,'String',0.4);

set(handles.corr_lag,'String',50);
set(handles.delay,'String',0);

% initialize the start path for browsing
handles.path_browse_root = '';

%fill the post proc option drop down list
vars{1} = 'Statistical time series analysis';
vars{2} = 'Generate activity maps';
vars{3} = 'Correlation analysis'; 
vars{4} = 'Scatter plots'; 
vars{5} = 'Protrusion vs. Retraction statistics'; 
%vars{6} = 'FFT analysis'; 
%vars{7} = 'Power spectrum';
set(handles.post_proc_operations,'String',vars);

set(handles.trend   ,'String', {'Filtered Data','Noise part'});

  
set(handles.filter_method   ,'String', {'No filtering',...
                                         'Filter by Diff. (detrend)',...
                                         'Filter with 4 poly',...
                                         '5 polynom',...
                                         'Moving average',...
                                         'Gauss filtering',...
                                         'Smoothing spline filtering'});

project_array = cell(1,1);
project_array(1) = cellstr('-- Choose a project --');
% set the limited entries in the popup menu
project_array(2) = cellstr('Browse ..');
set(handles.project_array   ,'String', project_array);
handles.nr_project_array = 2;
handles.project_def_file = 0;

% file the other dropdown menues
prot_dir_array = cell(1,1);
prot_dir_array(1) = cellstr('-- Choose protrusion directory --');
prot_dir_array(2) = cellstr('Browse ..');
set(handles.prot_dir_array   ,'String', prot_dir_array);
handles.nr_prot_dir_array = 2;

flow_dir_array = cell(1,1);
flow_dir_array(1) = cellstr('-- Choose mpm or vector file --');
flow_dir_array(2) = cellstr('Browse ..');
set(handles.flow_dir_array   ,'String', flow_dir_array);
handles.nr_flow_dir_array = 2;

scores_dir_array = cell(1,1);
scores_dir_array(1) = cellstr('-- Choose scores directory --');
scores_dir_array(2) = cellstr('Browse ..');
set(handles.scores_dir_array   ,'String', scores_dir_array);
handles.nr_scores_dir_array = 2;

activity_dir_array1 = cell(1,1);
activity_dir_array1(1) = cellstr('-- Choose image directory --');
activity_dir_array1(2) = cellstr('Browse ..');
set(handles.activity_dir_array1   ,'String', activity_dir_array1);
handles.nr_activity_dir_array1 = 2;

activity_dir_array2 = cell(1,1);
activity_dir_array2(1) = cellstr('-- Choose image directory --');
activity_dir_array2(2) = cellstr('Browse ..');
set(handles.activity_dir_array2   ,'String', activity_dir_array2);
handles.nr_activity_dir_array2 = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Post processing dropdown  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
merg_dir_array = cell(1,1);
merg_dir_array(1) = cellstr('Browse ..');
set(handles.merg_dir_array   ,'String', merg_dir_array);
handles.nr_merg_dir_array = 1;

var1_file_array = cell(1,1);
var1_file_array(1) = cellstr('Browse ..');
set(handles.var1_file_array   ,'String', var1_file_array);
handles.nr_var1_file_array = 1;

var2_file_array = cell(1,1);
var2_file_array(1) = cellstr('Browse ..');
set(handles.var2_file_array   ,'String', var2_file_array);
handles.nr_var2_file_array = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.pixel_per_frame(1),'Check','off');
set(handles.nm_per_sec(1),'Check','on');
set(handles.um_per_min(1),'Check','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.sum_scores(1),'Check','off');
set(handles.sum_scores_sec(1),'Check','on');
set(handles.av_scores(1),'Check','off');
set(handles.av_scores_sec(1),'Check','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.shift_scores_half_time_step(1),'Check','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.scores_pos(1),'Enable','off');
set(handles.scores_neg(1),'Enable','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disable the post processing functionalities
set(handles.post_proc_operations,'Enable','off');
set(handles.pralpha_postproc,'Enable','off');
set(handles.var1_file_array,'Enable','off');
set(handles.var2_file_array,'Enable','off');
set(handles.trend   ,'Enable', 'off');
set(handles.windowSize,'Enable','off');
set(handles.spline_s,'Enable','off');
set(handles.filter_method,'Enable','off');
set(handles.corr_lag,'Enable','off');
set(handles.delay,'Enable','off');
handles.interpolation = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default parameters for 'Time window and alignment'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.defProtTimeStep     = 1;
handles.defProtTimeWinL     = 3;
handles.defFlowTimeStep     = 1;
handles.defFlowTimeWinL     = 3;
handles.defScoreTimeStep    = 1;
handles.defScoreTimeWinL    = 1;
handles.defActivityTimeStep = 1;
handles.defActivityTimeWinL = 1;

handles.protTimeStep     = handles.defProtTimeStep;
handles.protTimeWinL     = handles.defProtTimeWinL;
handles.flowTimeStep     = handles.defFlowTimeStep;
handles.flowTimeWinL     = handles.defFlowTimeWinL;
handles.scoreTimeStep    = handles.defScoreTimeStep;
handles.scoreTimeWinL    = handles.defScoreTimeWinL;
handles.activityTimeStep = handles.defActivityTimeStep;
handles.activityTimeWinL = handles.defActivityTimeWinL;

handles.segShiftDir = 'normal';

% Choose default command line output for prMergePanel
handles.output = hObject;
guidata(hObject, handles);




function varargout = prMergePanel_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;



function total_time_steps_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function total_time_steps_Callback(hObject, eventdata, handles)


function first_time_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function first_time_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%      Merging       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_pr_alpha_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%     Determine which jobs have to be done   %%%%%%%%%%%%%%%%%%%%%
do_prot          = get(handles.do_prot,'Value');
do_flow          = get(handles.do_flow,'Value');
do_scores        = get(handles.do_scores,'Value');
do_activity_1    = get(handles.do_activity_1,'Value');
do_activity_2    = get(handles.do_activity_2,'Value');

%      fill that inromatiopn into the merge_paramters strucutre   %%%%%%%%%
merge_parameters.do_prot        = do_prot;
merge_parameters.do_flow        = do_flow;
merge_parameters.do_scores      = do_scores;
merge_parameters.do_activity_1  = do_activity_1;
merge_parameters.do_activity_2  = do_activity_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%    Determine the directories needed to run merge %%%%%%%%%%
% Get the project directory
val                             = get(handles.project_array,'Value');
project_dir_list                = get(handles.project_array,'String');
merge_parameters.project_dir    = project_dir_list{val};

% Get the protrusion directory
val             = get(handles.prot_dir_array,'Value');
prot_dir_list   = get(handles.prot_dir_array,'String');
prot_dir        = prot_dir_list{val};
merge_parameters.prot_dir = prot_dir; 


% Get the flow file
if do_flow
    val             = get(handles.flow_dir_array,'Value');
    flow_dir_list   = get(handles.flow_dir_array,'String');
    flow_dir        = flow_dir_list{val};

    merge_parameters.flow_dir = flow_dir;
else
    merge_parameters.flow_dir = '0';
end

% Get the scores file
if do_scores
    val             = get(handles.scores_dir_array,'Value');
    scores_dir_list = get(handles.scores_dir_array,'String');
    scores_dir      = scores_dir_list{val};

    merge_parameters.scores_dir = scores_dir;
else
    merge_parameters.scores_dir = '0';
end

% Get the activity 1 directory
if do_activity_1
    val                  = get(handles.activity_dir_array1,'Value');
    activity_dir_list    = get(handles.activity_dir_array1,'String');
    activity_dir         = activity_dir_list{val};

    merge_parameters.activity_dir_1 = activity_dir;
else
    merge_parameters.activity_dir_1 = '0';
end

% Get the activity 2 directory
if do_activity_2
    val                  = get(handles.activity_dir_array2,'Value');
    activity_dir_list    = get(handles.activity_dir_array2,'String');
    activity_dir         = activity_dir_list{val};

    merge_parameters.activity_dir_2 = activity_dir;
else
    merge_parameters.activity_dir_2 = '0';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check if to shift the scores and activity values by half
% time step to align it with protrusion measurements
% That is the default for assembly and disassembly scores.
if strcmp(get(handles.shift_scores_half_time_step(1),'Check'),'on');
    merge_parameters.shift_scores_half_time_step = 1;
else
    merge_parameters.shift_scores_half_time_step = 0;
end

% check what scores values to use
if strcmp(get(handles.scores_normal(1),'Check'),'on');
    merge_parameters.scores_sign = 0;
elseif strcmp(get(handles.scores_pos(1),'Check'),'on');
    merge_parameters.scores_sign = 1;
else
    merge_parameters.scores_sign = 2;
end

% check how to convert score values
if strcmp(get(handles.sum_scores(1),'Check'),'on');
    merge_parameters.scores_convert = 0;
elseif strcmp(get(handles.sum_scores_sec(1),'Check'),'on');
    merge_parameters.scores_convert = 1;
elseif strcmp(get(handles.av_scores(1),'Check'),'on');
    merge_parameters.scores_convert = 2;   
elseif strcmp(get(handles.av_scores_sec(1),'Check'),'on');
    merge_parameters.scores_convert = 3;
end

% get the data format of the score and flow values
if strcmp(get(handles.scores_matrix(1),'Check'),'on');
    merge_parameters.scores_matrix = 1;
else
    merge_parameters.scores_matrix = 0;
end

if strcmp(get(handles.flow_matrix(1),'Check'),'on');
    merge_parameters.flow_matrix = 1;
else
    merge_parameters.flow_matrix = 0;
end

% Get what unit system to use
if strcmp(get(handles.pixel_per_frame(1),'Check'),'on');
    merge_parameters.units = 0;
elseif strcmp(get(handles.nm_per_sec(1),'Check'),'on');
    merge_parameters.units = 1;    
elseif strcmp(get(handles.um_per_min(1),'Check'),'on');
    merge_parameters.units = 2;    
end

% Check if it is "plot vector" run
merge_parameters.plot_vectors = get(handles.plot_vectors, 'Value');

% Segmentation parameters   
seg_def = get(handles.seg_def,'Value');
if seg_def == 1
    merge_parameters.seg_nr     = str2num(get(handles.seg_val, 'String'));
    merge_parameters.seg_length = 0;
else
    merge_parameters.seg_nr     = 0;
    merge_parameters.seg_length = str2num(get(handles.seg_val, 'String'));
end

if strcmp(get(handles.dynamic_windows(1),'Check'),'on');
    merge_parameters.windows_type = 1;
elseif strcmp(get(handles.static_window(1),'Check'),'on');
    merge_parameters.windows_type = 2;    
elseif strcmp(get(handles.static_window_band(1),'Check'),'on');
    merge_parameters.windows_type = 3;
else
    merge_parameters.windows_type = 1;  
end

% Read the protrusion parameter file
prot_parameters = read_parameters([merge_parameters.project_dir filesep merge_parameters.prot_dir filesep 'parameters.dat'], 0);
if prot_parameters.contr == -99
    msgbox('Could not find any protrusion parameter file','Warning','warn');
    return
end
handles.prot_parameters = prot_parameters;

% get the image file name
[pathstr,name,ext,versn] = fileparts(prot_parameters.file);
% Check if LINUX or WIN system was used to generate edges
msgstr1 = 'Edge was run on a different system than the current system. ';
msgstr2 = 'Check the parameters.dat file in the edge directory. If the image ';
msgstr3 = 'path has forward slashes it was a LINUX otherwise a WIN system. ';
msgstr4 = 'You can: 1) run the merger on the correct system, or 2) manually ';
msgstr5 = 'edit the project and image path';
if ~isempty(findstr(name,':'));
    msgbox([msgstr1, msgstr2, msgstr3, msgstr4, msgstr5], 'System incompability','warn');
    return;
end

merge_parameters.img_name = name;

merge_parameters.first_time         = str2num(get(handles.first_time, 'String'));
merge_parameters.total_time_steps   = str2num(get(handles.total_time_steps, 'String'));
merge_parameters.start_seg          = str2num(get(handles.start_seg, 'String'));
merge_parameters.end_seg            = str2num(get(handles.end_seg, 'String'));
merge_parameters.seg_shift          = str2num(get(handles.seg_shift, 'String'));
merge_parameters.interpolation      = handles.interpolation;
merge_parameters.seg_depth          = str2num(get(handles.seg_depth, 'String'));
merge_parameters.d0                 = str2num(get(handles.d0, 'String'));
merge_parameters.replace_nan        = get(handles.replace_nan, 'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pass parameters for time window and alignment.
merge_parameters.protTimeStep = handles.protTimeStep;
merge_parameters.protTimeWinL = handles.protTimeWinL;
merge_parameters.flowTimeStep       = handles.flowTimeStep;
merge_parameters.flowTimeWinL       = handles.flowTimeWinL;
merge_parameters.scoreTimeStep      = handles.scoreTimeStep;
merge_parameters.scoreTimeWinL      = handles.scoreTimeWinL;
merge_parameters.activityTimeStep   = handles.activityTimeStep;
merge_parameters.activityTimeWinL   = handles.activityTimeWinL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

merge_parameters.segShiftDir   = handles.segShiftDir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the prMerge function   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prMerge(merge_parameters, prot_parameters);

guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     Post processing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pralpha_postproc_Callback(hObject, eventdata, handles)

post_parameters.post_proc_options = get(handles.post_proc_operations,'Value');

% get the project directory
val                         = get(handles.project_array,'Value');
project_dir_list            = get(handles.project_array,'String');
post_parameters.project_dir = project_dir_list{val};

% get the merge directory
val                         = get(handles.merg_dir_array,'Value');
post_parameters.merg_dir    = handles.merg_dir{val};

 % get the protrusion directory
val                         = get(handles.prot_dir_array,'Value');
prot_dir_list               = get(handles.prot_dir_array,'String');
post_parameters.prot_dir    = prot_dir_list{val}; 

% get the var 1 file
val                         = get(handles.var1_file_array,'Value');
post_parameters.var1_file   = handles.var1_file{val};

% get the var 2 file
val                         = get(handles.var2_file_array,'Value');
post_parameters.var2_file   = handles.var2_file{val};

% check if data should be filtered
val                       = get(handles.filter_method,'Value');
post_parameters.do_filter = val;

% get the window size of the filter
post_parameters.windowSize = str2num(get(handles.windowSize,'String'));
post_parameters.spline_s   = str2num(get(handles.spline_s,'String'));

% check if the filtered data should be processed or the noise
% content
val                   = get(handles.trend,'Value');
post_parameters.trend = val;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Operation specific data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the lag for correlation
post_parameters.corr_lag = str2num(get(handles.corr_lag,'String'));

% get the delay for the variable space plot
post_parameters.delay = str2num(get(handles.delay,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%  save the post_parameters in merge directory %%%%%%%%%%%%%%%%%%%%%%%%
post_parFile = [post_parameters.merg_dir filesep 'post_parameters.mat'];
save(post_parFile,'post_parameters');

% run the post processing
prMergePostProc(post_parameters);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function protrusion_Callback(hObject, eventdata, handles)
prot_value = get(hObject,'Value');

if prot_value
    set(handles.write_prot_data,'Enable','on');
    set(handles.write_prot_data,'Value',1);   
else
    set(handles.write_prot_data,'Value',0);
    set(handles.write_prot_data,'Enable','off');
end
guidata(hObject, handles);


function plot_vectors_Callback(hObject, eventdata, handles)

function file_Callback(hObject, eventdata, handles)
function exit_Callback(hObject, eventdata, handles)
delete(handles.figure1);

function help_Callback(hObject, eventdata, handles)

function about_Callback(hObject, eventdata, handles)

function about_menu_Callback(hObject, eventdata, handles)
prAboutMerge;
function help_menu_Callback(hObject, eventdata, handles)
function helph_menu_Callback(hObject, eventdata, handles)
dir_path_full = which('prPanel');
[pathstr,name,ext,versn] = fileparts(dir_path_full);
web (['file:///', pathstr, filesep, 'prHelpMain.html']); 


function project_array_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Project callback   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function project_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
merge_parameters.project_dir = char(string_list{val}); 


if strcmp(merge_parameters.project_dir,'Browse ..')
    dir_w = uigetdir(handles.path_browse_root,'Select/create directory for results');
    if dir_w == 0
        return;
    end
    handles.path_browse_root = dir_w;
    string_list(handles.nr_project_array)   = cellstr(dir_w);
    string_list(handles.nr_project_array+1) = cellstr('Browse ..');
    
    
    set(handles.project_array   ,'String', string_list);
    handles.nr_project_array = handles.nr_project_array+1;
    guidata(hObject, handles);  
    
    merge_parameters.project_dir = dir_w;

elseif handles.project_def_file == 0
   val = get(hObject,'Value');
   string_list = get(hObject,'String');
   dir_w = char(string_list{val});   
end

% search the directory for edge sub directories
prot_dir = dir([dir_w filesep 'edge*']);
% fill the menu
nr_prot_dir = size(prot_dir,1);
if nr_prot_dir > 0
    set(handles.prot_dir_array,'Value',nr_prot_dir);
    for i = 1:nr_prot_dir
        prot_dir_array(i) = cellstr(prot_dir(i).name);
    end
    set(handles.prot_dir_array ,'String', prot_dir_array);
else
    set(handles.prot_dir_array,'Value',1);
    set(handles.prot_dir_array,'String', 'No directory found');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prot_dir_array_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function prot_dir_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
merge_parameters.prot_dir = char(string_list{val}); 


if strcmp(merge_parameters.prot_dir,'Browse ..')
    dir_w = uigetdir('start_path','Select directory containing protrusion');
    if dir_w == 0
        return;
    end
    string_list(handles.nr_prot_dir_array)   = cellstr(dir_w);
    string_list(handles.nr_prot_dir_array+1) = cellstr('Browse ..');

    set(handles.prot_dir_array   ,'String', string_list);
    handles.nr_prot_dir_array = handles.nr_prot_dir_array+1;
    guidata(hObject, handles);

    merge_parameters.prot_dir = dir_w;
end

function do_prot_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flow_dir_array_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function flow_dir_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
merge_parameters.flow_dir = char(string_list{val}); 


if strcmp(merge_parameters.flow_dir,'Browse ..')
    [file_n dir_w_n] = uigetfile(handles.path_browse_root,'Select mpm file or first vector file');
    dir_w = [dir_w_n file_n];
    if dir_w == 0
        return;
    end
    string_list(handles.nr_flow_dir_array)   = cellstr(dir_w);
    string_list(handles.nr_flow_dir_array+1) = cellstr('Browse ..');

    set(handles.flow_dir_array   ,'String', string_list);
    handles.nr_flow_dir_array = handles.nr_flow_dir_array+1;
    guidata(hObject, handles);

    merge_parameters.flow_dir = dir_w;
end

function do_flow_Callback(hObject, eventdata, handles)
ans = get(handles.do_flow,'Value');
if ans == 1
    set(handles.d0,'Enable','on');
else
    set(handles.d0,'Enable','off');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scores_dir_array_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function scores_dir_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
merge_parameters.scores_dir = char(string_list{val}); 

if strcmp(merge_parameters.scores_dir,'Browse ..')
    [file_n dir_w_n] = uigetfile(handles.path_browse_root,'Select file containing scores');
    dir_w = [dir_w_n file_n];
    if dir_w == 0
        return;
    end
    string_list(handles.nr_scores_dir_array)   = cellstr(dir_w);
    string_list(handles.nr_scores_dir_array+1) = cellstr('Browse ..');

    set(handles.scores_dir_array   ,'String', string_list);
    handles.nr_scores_dir_array = handles.nr_scores_dir_array+1;
    guidata(hObject, handles);

    merge_parameters.scores_dir = dir_w;
end

function do_scores_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function activity_dir_array1_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
merge_parameters.activity1_dir = char(string_list{val}); 

if strcmp(merge_parameters.activity1_dir,'Browse ..')
    dir_w = uigetdir(handles.path_browse_root,'Select directory with images');
    if dir_w == 0
        return;
    end
    string_list(handles.nr_activity_dir_array1)   = cellstr(dir_w);
    string_list(handles.nr_activity_dir_array1+1) = cellstr('Browse ..');

    set(handles.activity_dir_array1   ,'String', string_list);
    handles.nr_activity_dir_array1 = handles.nr_activity_dir_array1+1;
    guidata(hObject, handles);

    merge_parameters.activity1_dir = dir_w;
end


function activity_dir_array1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function activity_dir_array2_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
merge_parameters.activity2_dir = char(string_list{val}); 

if strcmp(merge_parameters.activity2_dir,'Browse ..')
    dir_w = uigetdir(handles.path_browse_root,'Select directory with images');
    if dir_w == 0
        return;
    end
    string_list(handles.nr_activity_dir_array2)   = cellstr(dir_w);
    string_list(handles.nr_activity_dir_array2+1) = cellstr('Browse ..');

    set(handles.activity_dir_array2   ,'String', string_list);
    handles.nr_activity_dir_array2 = handles.nr_activity_dir_array2+1;
    guidata(hObject, handles);

    merge_parameters.activity2_dir = dir_w;
end


function activity_dir_array2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function do_activity_1_Callback(hObject, eventdata, handles)
function do_activity_2_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function start_seg_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function start_seg_Callback(hObject, eventdata, handles)


function end_seg_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function end_seg_Callback(hObject, eventdata, handles)


function seg_shift_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function seg_shift_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Determine which post processing operation to run  %%%%%%%

function post_proc_operations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function post_proc_operations_Callback(hObject, eventdata, handles)
% vars{1} = 'Statistical time series analysis';
% vars{2} = 'Generate activity maps';
% vars{3} = 'Correlation analysis'; 
% vars{4} = 'Scatter plots'; 
% vars{5} = 'Protrusion vs. Retraction statistics'; 
% vars{6} = 'FFT analysis'; 
% vars{7} = 'Power spectrum';
val = get(hObject,'Value');
if val == 3
    set(handles.corr_lag,'Enable','on');
elseif val == 4 || val == 5
    set(handles.delay,'Enable','on');
else
    set(handles.delay,'Enable','off'); 
    set(handles.corr_lag,'Enable','off');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function merg_dir_array_Callback(hObject, eventdata, handles)
% get the project directory
val                 = get(handles.project_array,'Value');
project_dir_list    = get(handles.project_array,'String');
project_dir         = project_dir_list{val};

val              = get(handles.merg_dir_array,'Value');
merg_dir_list    = get(handles.merg_dir_array,'String');
merg_dir         = merg_dir_list{val};

handles.variable_path = [project_dir filesep merg_dir];
guidata(hObject, handles);





function merg_dir_array_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function var1_file_array_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function var1_file_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
browse_ans = char(string_list{val}); 

if strcmp(browse_ans,'Browse ..')
    if isfield(handles,'variable_path')
        cd(handles.variable_path);
    end
    [name, path] = uigetfile('*.mat','Select merged data set file');
    if path == 0
        return;
    end 
    filsep_pos =  strfind(path, filesep);
    part_path = path(filsep_pos(end-2):end);
    
    part_path = [part_path name];
    string_list(handles.nr_var1_file_array)   = cellstr(part_path);
    string_list(handles.nr_var1_file_array+1) = cellstr('Browse ..');

    handles.var1_file{handles.nr_var1_file_array} = [path name];
    
    set(handles.var1_file_array,'String', string_list);
    handles.nr_var1_file_array = handles.nr_var1_file_array+1;
    guidata(hObject, handles); 
end

function var2_file_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
browse_ans = char(string_list{val}); 

if strcmp(browse_ans,'Browse ..')
    if isfield(handles,'variable_path')
        cd(handles.variable_path);
    end
    [name, path] = uigetfile('*.mat','Select merged data set file');
    if path == 0
        return;
    end 
    filsep_pos =  strfind(path, filesep);
    part_path = path(filsep_pos(end-2):end);
    
    part_path = [part_path name];   
    string_list(handles.nr_var2_file_array)   = cellstr(part_path);
    string_list(handles.nr_var2_file_array+1) = cellstr('Browse ..');

    handles.var2_file{handles.nr_var2_file_array} = [path name];
    
    set(handles.var2_file_array,'String', string_list);
    handles.nr_var2_file_array = handles.nr_var2_file_array+1;
    guidata(hObject, handles);
end



function var2_file_array_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function pr_alpha_option_Callback(hObject, eventdata, handles)
set(handles.pr_alpha_option(1),'Check','on');
set(handles.post_proc_option(1),'Check','off');

% disable the post processing functionalities
set(handles.post_proc_operations,'Enable','off');
set(handles.pralpha_postproc,'Enable','off');
set(handles.var1_file_array,'Enable','off');
set(handles.var2_file_array,'Enable','off');
set(handles.filter_method,'Enable','off');
set(handles.trend   ,'Enable', 'off');
set(handles.windowSize,'Enable','off');
set(handles.spline_s,'Enable','off');
set(handles.corr_lag,'Enable','off');
set(handles.delay,'Enable','off');

% enable the Merge functionalities
set(handles.first_time,'Enable','on');
set(handles.total_time_steps,'Enable','on');
set(handles.start_seg,'Enable','on');
set(handles.end_seg,'Enable','on');
set(handles.seg_shift,'Enable','on');
set(handles.run_pr_alpha,'Enable','on');
set(handles.plot_vectors,'Enable','on');
if get(handles.do_flow,'Value');
    set(handles.d0,'Enable','on');
end

set(handles.seg_val,'Enable','on');
set(handles.seg_def,'Enable','off');
set(handles.seg_depth,'Enable','on');

% --------------------------------------------------------------------
function post_proc_option_Callback(hObject, eventdata, handles)
set(handles.pr_alpha_option(1),'Check','off');
set(handles.post_proc_option(1),'Check','on');
% enable the post processing functionalities
set(handles.post_proc_operations,'Enable','on');
set(handles.pralpha_postproc,'Enable','on');
set(handles.var1_file_array,'Enable','on');
set(handles.var2_file_array,'Enable','on');

set(handles.filter_method,'Enable','on');
set(handles.filter_method,'Value',1);
set(handles.trend   ,'Enable', 'off');
set(handles.windowSize,'Enable','off');
set(handles.spline_s,'Enable','off');
set(handles.corr_lag,'Enable','on');
set(handles.delay,'Enable','on');

% get the project directory
val                 = get(handles.project_array,'Value');
project_dir_list    = get(handles.project_array,'String');
dir_w               = project_dir_list{val};

% Fill drop down menu

% Search the directory for merg projects
merg_dir = dir([dir_w filesep 'merg*']);
% fill the menu
nr_merg_dir = size(merg_dir,1);
if nr_merg_dir > 0
    %set(handles.merg_dir_array,'Value',nr_merg_dir);
    for i = 1:nr_merg_dir
        merg_dir_name(i) = cellstr(merg_dir(i).name);
        tmp = [dir_w filesep char(merg_dir_name(i))];
        handles.merg_dir{i} = tmp;
    end
    nr_el = length(get(handles.merg_dir_array ,'String'));
    for i = nr_merg_dir+1 : nr_el
        handles.merg_dir{i} = '-'; 
    end
    set(handles.merg_dir_array ,'String', merg_dir_name);
    %string_list(handles.nr_var1_file_array+1) = cellstr('Browse ..');
end

% disable the prMerge functionalities
set(handles.first_time,'Enable','off');
set(handles.total_time_steps,'Enable','off');
set(handles.start_seg,'Enable','off');
set(handles.end_seg,'Enable','off');
set(handles.seg_shift,'Enable','off');
set(handles.run_pr_alpha,'Enable','off');

set(handles.seg_val,'Enable','off');
set(handles.seg_def,'Enable','off');
set(handles.seg_depth,'Enable','off');

set(handles.plot_vectors,'Enable','off');

set(handles.d0,'Enable','off');
guidata(hObject, handles);



function no_interpolation_Callback(hObject, eventdata, handles)
handles.interpolation = 0;
set(handles.no_interpolation(1),'Check','on');
set(handles.dir_interpolation(1),'Check','off');
set(handles.s4_interpolation(1),'Check','off');
guidata(hObject, handles);

function dir_interpolation_Callback(hObject, eventdata, handles)
handles.interpolation = 2;
set(handles.no_interpolation(1),'Check','off');
set(handles.dir_interpolation(1),'Check','on');
set(handles.s4_interpolation(1),'Check','off');
guidata(hObject, handles);

function s4_interpolation_Callback(hObject, eventdata, handles)
handles.interpolation = 1;
set(handles.no_interpolation(1),'Check','off');
set(handles.dir_interpolation(1),'Check','off');
set(handles.s4_interpolation(1),'Check','on');
guidata(hObject, handles);


function filter_method_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if val==1
    set(handles.trend,'Enable','off');
    set(handles.spline_s,'Enable','off');
    set(handles.windowSize,'Enable','off');
elseif val == 2 | val == 3 | val == 4
    set(handles.trend,'Enable','on');
    set(handles.spline_s,'Enable','off');
    set(handles.windowSize,'Enable','off');
elseif val == 5 | val == 6
    set(handles.trend,'Enable','on');
    set(handles.spline_s,'Enable','off');    
    set(handles.windowSize,'Enable','on');
elseif val == 7
    set(handles.trend,'Enable','on');
    set(handles.spline_s,'Enable','on');    
    set(handles.windowSize,'Enable','off');      
end


function filter_method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowSize_Callback(hObject, eventdata, handles)
function windowSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function trend_Callback(hObject, eventdata, handles)
function trend_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function spline_s_Callback(hObject, eventdata, handles)
function spline_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function corr_lag_Callback(hObject, eventdata, handles)
function corr_lag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function seg_nr_Callback(hObject, eventdata, handles)
function seg_nr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function seg_depth_Callback(hObject, eventdata, handles)
function seg_depth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function seg_val_Callback(hObject, eventdata, handles)
function seg_val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function seg_def_Callback(hObject, eventdata, handles)
function seg_def_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function delay_Callback(hObject, eventdata, handles)
function delay_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_vector_s_Callback(hObject, eventdata, handles)
set(handles.plot_vector_s,'Value',1);
set(handles.plot_vector_m,'Value',0);

function plot_vector_m_Callback(hObject, eventdata, handles)
set(handles.plot_vector_s,'Value',0);
set(handles.plot_vector_m,'Value',1);


function colormap_Callback(hObject, eventdata, handles)
handles.color_map_parameters = prColormapPanel;
guidata(hObject, handles);

function parameters_Callback(hObject, eventdata, handles)



function d0_Callback(hObject, eventdata, handles)

function d0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scores_settings_Callback(hObject, eventdata, handles)

function scores_normal_Callback(hObject, eventdata, handles)
set(handles.scores_normal(1),'Check','on');
set(handles.scores_pos(1),'Check','off');
set(handles.scores_neg(1),'Check','off');

function scores_pos_Callback(hObject, eventdata, handles)
set(handles.scores_normal(1),'Check','off');
set(handles.scores_pos(1),'Check','on');
set(handles.scores_neg(1),'Check','off');

function scores_neg_Callback(hObject, eventdata, handles)
set(handles.scores_normal(1),'Check','off');
set(handles.scores_pos(1),'Check','off');
set(handles.scores_neg(1),'Check','on');

function scores_data_format_Callback(hObject, eventdata, handles)
function scores_matrix_Callback(hObject, eventdata, handles)
set(handles.scores_matrix(1),'Check','on');
set(handles.scores_vector(1),'Check','off');

function scores_vector_Callback(hObject, eventdata, handles)
set(handles.scores_matrix(1),'Check','off');
set(handles.scores_vector(1),'Check','on');

function flow_data_format_Callback(hObject, eventdata, handles)
function flow_matrix_Callback(hObject, eventdata, handles)
set(handles.flow_matrix(1),'Check','on');
set(handles.flow_vector(1),'Check','off');

function flow_vector_Callback(hObject, eventdata, handles)
set(handles.flow_matrix(1),'Check','off');
set(handles.flow_vector(1),'Check','on');

function interpolation_Callback(hObject, eventdata, handles)

function units_Callback(hObject, eventdata, handles)

function pixel_per_frame_Callback(hObject, eventdata, handles)
set(handles.pixel_per_frame(1),'Check','on');
set(handles.nm_per_sec(1),'Check','off');
set(handles.um_per_min(1),'Check','off');

function nm_per_sec_Callback(hObject, eventdata, handles)
set(handles.pixel_per_frame(1),'Check','off');
set(handles.nm_per_sec(1),'Check','on');
set(handles.um_per_min(1),'Check','off');

function um_per_min_Callback(hObject, eventdata, handles)
set(handles.pixel_per_frame(1),'Check','off');
set(handles.nm_per_sec(1),'Check','off');
set(handles.um_per_min(1),'Check','on');

function replace_nan_Callback(hObject, eventdata, handles)

function score_av_Callback(hObject, eventdata, handles)

function sum_scores_Callback(hObject, eventdata, handles)
set(handles.sum_scores(1),'Check','on');
set(handles.av_scores(1),'Check','off');
set(handles.av_scores_sec(1),'Check','off');
set(handles.sum_scores_sec(1),'Check','off');

function av_scores_Callback(hObject, eventdata, handles)
set(handles.sum_scores(1),'Check','off');
set(handles.av_scores(1),'Check','on');
set(handles.av_scores_sec(1),'Check','off');
set(handles.sum_scores_sec(1),'Check','off');

function av_scores_sec_Callback(hObject, eventdata, handles)
set(handles.sum_scores(1),'Check','off');
set(handles.av_scores(1),'Check','off');
set(handles.av_scores_sec(1),'Check','on');
set(handles.sum_scores_sec(1),'Check','off');

function sum_scores_sec_Callback(hObject, eventdata, handles)
set(handles.sum_scores(1),'Check','off');
set(handles.av_scores(1),'Check','off');
set(handles.av_scores_sec(1),'Check','off');
set(handles.sum_scores_sec(1),'Check','on');

function shift_scores_half_time_step_Callback(hObject, eventdata, handles)
if strcmp(get(handles.shift_scores_half_time_step(1),'Check'),'on');
    set(handles.shift_scores_half_time_step(1),'Check','off');
else
    set(handles.shift_scores_half_time_step(1),'Check','on');
end




function dynamic_windows_Callback(hObject, eventdata, handles)
set(handles.dynamic_windows(1),'Check','on');
set(handles.static_window(1),'Check','off');
set(handles.static_window_band(1),'Check','off');

function static_window_Callback(hObject, eventdata, handles)
set(handles.dynamic_windows(1),'Check','off');
set(handles.static_window(1),'Check','on');
set(handles.static_window_band(1),'Check','off');

function static_window_band_Callback(hObject, eventdata, handles)
set(handles.dynamic_windows(1),'Check','off');
set(handles.static_window(1),'Check','off');
set(handles.static_window_band(1),'Check','on');




% --------------------------------------------------------------------
function timeAlign_parameter_Callback(hObject, eventdata, handles)
% hObject    handle to timeAlign_parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


defProtTimeStep  = handles.defProtTimeStep;
defProtTimeWinL  = handles.defProtTimeWinL;
defFlowTimeStep        = handles.defFlowTimeStep;
defFlowTimeWinL        = handles.defFlowTimeWinL;
defScoreTimeStep       = handles.defScoreTimeStep;
defScoreTimeWinL       = handles.defScoreTimeWinL;
defActivityTimeStep    = handles.defActivityTimeStep;
defActivityTimeWinL    = handles.defActivityTimeWinL;


protTimeStepStr = sprintf(['Unit: frame.\n\n' ...
   'Protrusion setting:\n' ...
   '  Enter the time step (default: %d)'], defProtTimeStep);
protTimeWinLStr = sprintf('  Enter the time window length (default: %d)', ...
   defProtTimeWinL);
flowTimeStepStr = sprintf(['Flow setting:\n' ...
   '  Enter the time step (default: %d)'], defFlowTimeStep);
flowTimeWinLStr = sprintf('  Enter the time window length (default: %d)', ...
   defFlowTimeWinL);
scoreTimeStepStr = sprintf(['Score setting:\n' ...
   '  Enter the time step (default: %d)'], defScoreTimeStep);
scoreTimeWinLStr = sprintf('  Enter the time window length (default: %d)', ...
   defScoreTimeWinL);
actTimeStepStr = sprintf(['Activity map setting:\n' ...
   '  Enter the time step (default: %d)'], defActivityTimeStep);
actTimeWinLStr = sprintf('  Enter the time window length (default: %d)', ...
   defActivityTimeWinL);
ans = inputdlg({protTimeStepStr, protTimeWinLStr, ...
   flowTimeStepStr, flowTimeWinLStr, ...
   scoreTimeStepStr, scoreTimeWinLStr, ...
   actTimeStepStr, actTimeWinLStr},'Parameters of Time Alignment',1, ...
   {num2str(defProtTimeStep),num2str(defProtTimeWinL), ...
   num2str(defFlowTimeStep),num2str(defFlowTimeWinL), ...
   num2str(defScoreTimeStep),num2str(defScoreTimeWinL), ...
   num2str(defActivityTimeStep),num2str(defActivityTimeWinL), ...
   });

if isempty(ans)
   return;
end

handles.protTimeStep     = str2num(ans{1});
handles.protTimeWinL     = str2num(ans{2});
handles.flowTimeStep     = str2num(ans{3});
handles.flowTimeWinL     = str2num(ans{4});
handles.scoreTimeStep    = str2num(ans{5});
handles.scoreTimeWinL    = str2num(ans{6});
handles.activityTimeStep = str2num(ans{7});
handles.activityTimeWinL = str2num(ans{8});

guidata(hObject, handles);




% --------------------------------------------------------------------
function normal_shift_Callback(hObject, eventdata, handles)
% hObject    handle to normal_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.normal_shift,'checked','on');
set(handles.parallel_shift,'checked','off');

handles.segShiftDir = 'normal';

guidata(hObject, handles);


% --------------------------------------------------------------------
function parallel_shift_Callback(hObject, eventdata, handles)
% hObject    handle to parallel_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.normal_shift,'checked','off');
set(handles.parallel_shift,'checked','on');
handles.segShiftDir = 'parallel';


guidata(hObject, handles);

