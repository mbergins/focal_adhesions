function varargout = prMergeProtrusionPanel(varargin)
% PRMERGEPROTRUSIONPANEL M-file for prMergeProtrusionPanel.fig
%      PRMERGEPROTRUSIONPANEL, by itself, creates a new PRMERGEPROTRUSIONPANEL or raises the existing
%      singleton*.
%
%      H = PRMERGEPROTRUSIONPANEL returns the handle to a new PRMERGEPROTRUSIONPANEL or the handle to
%      the existing singleton*.
%
%      PRMERGEPROTRUSIONPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRMERGEPROTRUSIONPANEL.M with the given input arguments.
%
%      PRMERGEPROTRUSIONPANEL('Property','Value',...) creates a new PRMERGEPROTRUSIONPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prMergeProtrusionPanel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prMergeProtrusionPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prMergeProtrusionPanel

% Last Modified by GUIDE v2.5 31-Jul-2006 15:29:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prMergeProtrusionPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @prMergeProtrusionPanel_OutputFcn, ...
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


% --- Executes just before prMergeProtrusionPanel is made visible.
function prMergeProtrusionPanel_OpeningFcn(hObject, eventdata, handles, varargin)

%fill the text boxes 
handles.parameter_set = 0;
set(handles.total_time_steps,'String',10);
handles.debug = 0;

% fill the segment definition list
set(handles.seg_val,'String',30);

% initialize the start path for browsing
handles.path_browse_root = '';


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.pixel_per_frame(1),'Check','off');
set(handles.nm_per_sec(1),'Check','on');
set(handles.um_per_min(1),'Check','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for prMergeProtrusionPanel
handles.output = hObject;
guidata(hObject, handles);




function varargout = prMergeProtrusionPanel_OutputFcn(hObject, eventdata, handles)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%    segmentation parameters   
merge_parameters.seg_nr     = str2num(get(handles.seg_val, 'String'));


% Read the protrusion parameter file
prot_parameter = read_parameters([merge_parameters.project_dir filesep merge_parameters.prot_dir filesep 'parameters.dat'], 0);
if prot_parameter.contr == -99
    msgbox('Could not find any protrusion parameter file','Warning','warn');
    return
end
handles.prot_parameter = prot_parameter;

% get the image file name
[pathstr,name,ext,versn] = fileparts(prot_parameter.file);
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


merge_parameters.total_time_steps   = str2num(get(handles.total_time_steps, 'String'));
merge_parameters.start_seg          = str2num(get(handles.start_seg, 'String'));
merge_parameters.end_seg            = str2num(get(handles.end_seg, 'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the prAlpha function   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prMergeProtrusion(merge_parameters, prot_parameter);
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

% search the directory for PROTRUSION projects
prot_dir = dir([dir_w filesep 'edge*']);
% fill the menue
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pr_alpha_dir_array_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function pr_alpha_dir_array_Callback(hObject, eventdata, handles)
% get the project directory
val                 = get(handles.project_array,'Value');
project_dir_list    = get(handles.project_array,'String');
project_dir         = project_dir_list{val};
% %get the pr alpha directory
val                  = get(handles.pr_alpha_dir_array,'Value');
pr_alpha_dir_list    = get(handles.pr_alpha_dir_array,'String');
pr_alpha_dir         = pr_alpha_dir_list{val};


% check for merge parameter file
fid = fopen([project_dir  pr_alpha_dir 'merge_parameters.dat'],'r');
if fid ~= -1
    ans = questdlg('Parameter file found. Take its values?');
    if strcmp(ans,'Yes')
        merge_parameters = read_pralpha_parameters([project_dir pr_alpha_dir 'merge_parameters.dat']);
        %write parameters to tmp disk
        %save_parameters([], parameters);
        handles.merge_parameters_set = 1; 
    end
    
    %set the GUI values accordingly
    set(handles.first_time,'String',num2str(merge_parameters.first_time));
    set(handles.total_time_steps,'String',num2str(merge_parameters.total_time_steps));
    set(handles.start_seg,'String',num2str(merge_parameters.start_seg));
    set(handles.end_seg,'String',num2str(merge_parameters.end_seg));
    set(handles.seg_shift,'String',num2str(merge_parameters.seg_shift));
end









function seg_nr_Callback(hObject, eventdata, handles)
function seg_nr_CreateFcn(hObject, eventdata, handles)
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

function parameters_Callback(hObject, eventdata, handles)



function scores_vector_Callback(hObject, eventdata, handles)
set(handles.scores_matrix(1),'Check','off');
set(handles.scores_vector(1),'Check','on');



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

function Untitled_4_Callback(hObject, eventdata, handles)

