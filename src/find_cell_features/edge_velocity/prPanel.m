function varargout = prPanel(varargin)
% PRPANEL M-file for prPanel.fig
%      PRPANEL, by itself, creates a new PRPANEL or raises the existing
%      singleton*.
%
%      H = PRPANEL returns the handle to a new PRPANEL or the handle to
%      the existing singleton*.
%
%      PRPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRPANEL.M with the given input arguments.
%
%      PRPANEL('Property','Value',...) creates a new PRPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prPanel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prPanel

% Last Modified by GUIDE v2.5 31-Jul-2006 11:05:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @prPanel_OutputFcn, ...
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


% --- Executes just before prPanel is made visible.
function prPanel_OpeningFcn(hObject, eventdata, handles, varargin)

%fill the text boxes 
handles.parameters_set = 0;
set(handles.max_img,'String',10);
set(handles.t_step,'String',1);
set(handles.first_image,'String',1);
parameters.debug = 0;
handles.parameters.debug = 0;

% initialize the start path for browsing
handles.path_browse_root = '';
handles.img_browse_root = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Check whether fsmCenter is running - if yes get project %%%%%%%%%%%
% from it
hfsmC=findall(0,'Tag','fsmCenter','Name','fsmCenter');
if ~isempty(hfsmC)
    % Get current project from fsmCenter
    handlesFsmCenter=guidata(hfsmC);
    settings=get(handlesFsmCenter.fsmCenter,'UserData');
    handles.fsmCenter_parameters.NA             = handlesFsmCenter.physiParam.NA;
    handles.fsmCenter_parameters.waveLength     = handlesFsmCenter.physiParam.waveLen;
    handles.fsmCenter_parameters.bitDepth       = handlesFsmCenter.physiParam.bitDepth;
    handles.fsmCenter_parameters.pixelSize      = handlesFsmCenter.physiParam.pixelSize; 
    handles.fsmCenter_parameters.frameInterval  = handlesFsmCenter.physiParam.frameInterval; 
    handles.fsmCenter_parameters.psfSigma       = handlesFsmCenter.physiParam.psfSigma; 
    
    set(handles.dir_array   ,'String', [settings.projDir filesep]);
    handles.fsmCenterEdgeDirectory = settings.subProjects(4);
    set(handles.image_array   ,'String', cellstr([char(settings.imageDir) char(settings.firstImgList(1))]));
else
    % set the initial options in the drop-down menu
    image_array = cell(1,1);
    image_array(1) = cellstr('-- Choose an image --');
    image_array(2) = cellstr('Browse ..');
    set(handles.image_array   ,'String', image_array);
    handles.nr_image_array = 2;

    %set the limited entries in the popup menu
    dir_array = cell(1,1);
    dir_array(1) = cellstr('-- Choose an working directory --');
    dir_array(2) = cellstr('Browse ..');
    set(handles.dir_array   ,'String', dir_array);
    handles.nr_dir_array = 2;
end

set(handles.img_depth   ,'String', {'8', '10', '12', '14', '16'});

% Choose default command line output for prPanel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


function varargout = prPanel_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function settings_Callback(hObject, eventdata, handles)
MAX_IMG             = str2num(get(handles.max_img,'String'));
T_STEP              = str2num(get(handles.t_step,'String'));
        

% if there is already a parameter file read this one
if handles.parameters_set == 1
    file_path = [];
    
    if ispc
        tmp_dir = getenv('TMP'); 
        tmp_dir = [tmp_dir filesep 'pr_parameters.mat'];        
    else
        tmp_dir = getenv('HOME');   
        tmp_dir = [tmp_dir filesep 'pr_parameters.mat'];        
    end
    % now load the parameters
    load(tmp_dir);
    
    % if the user changed to use bw images, re-set the filter
    % parameters
    if get(handles.use_bw_images,'Value');
        parameters.use_bw_images = 1;
        
        parameters.filter_image     = 0;
        parameters.erode_dilate     = 0;
        parameters.median_f         = 0;
        parameters.manual_thresh    = 1;
    else
        parameters.use_bw_images = 0;
    end
    
    
    % if fsmCenter is running take its parameters
    if isfield(handles,'fsmCenter_parameters')
        parameters.pixel                =   handles.fsmCenter_parameters.pixelSize;
        parameters.time_interval        =   handles.fsmCenter_parameters.frameInterval;
        fsmCenterRunning = 1;
    else
        fsmCenterRunning = 0;
    end      
    
    % call parameter panel with parameters to fill in
    [parameters] = prPanelPar('parameters',parameters,'status',fsmCenterRunning);
else
    % call the parameter panel default parameters
    
    % If we use already segmented images we have to use
    % specific settings
    if get(handles.use_bw_images,'Value');
        parameters.use_bw_images = 1;
        
        parameters.filter_image     = 0;
        parameters.erode_dilate     = 0;
        parameters.median_f         = 0;
        parameters.manual_thresh    = 1;
    else
        parameters.use_bw_images = 0;
        
        parameters.na               = 1.4;
        parameters.lambda           = 500;
        parameters.pixel            = 67;
        parameters.time_interval    = 10;
        parameters.prot_sampling    = 3;
        parameters.parenth_l        = 0;
        parameters.parenth_r        = 0;
        parameters.filter_image     = 1;
        parameters.img_sigma        = 0.9;
        parameters.erode_dilate     = 6;
        parameters.median_f         = 1;
        parameters.f_window         = 5;
        parameters.f_sigma          = 0.1;
        parameters.tolerance        = 20;
        parameters.mechanical       = 1;
        parameters.nearest          = 0;
        parameters.normal           = 0;
        parameters.tol              = 0;
        parameters.robust_min       = 1;
        parameters.k_s              = 0.01;
        parameters.k_w              = 1;
        parameters.cluster          = 0;
        parameters.cluster_method   = 'em';
        parameters.kmeans           = 0;
        parameters.em               = 1;
        parameters.manual_thresh    = 0;
        parameters.k_cluster        = 3;
        parameters.k_min            = 3;
        parameters.k_max            = 3;
        parameters.orient           = 0;
        parameters.cell_mode        = 1;
    end
    
    % if fsmCenter is running take its parameters
    if isfield(handles,'fsmCenter_parameters')
        parameters.pixel                =   handles.fsmCenter_parameters.pixelSize;
        parameters.time_interval        =   handles.fsmCenter_parameters.frameInterval;
        fsmCenterRunning = 1;
    else
        fsmCenterRunning = 0;
    end
     
    [parameters] = prPanelPar('parameters',parameters,'status',fsmCenterRunning);
end

%save parameters so they can be accessed by other callback functions
handles.parameters = parameters;
handles.parameters_set = 1;

% save the parameters so that they can be retrived 
% if the parameter panel is opnened again
if ispc
    tmp_dir = getenv('TMP');
    tmp_dir = [tmp_dir filesep 'pr_parameters.mat'];
else
    tmp_dir = getenv('HOME');
    tmp_dir = [tmp_dir filesep 'pr_parameters.mat'];
end
save(tmp_dir,'parameters');

guidata(hObject, handles);




function max_img_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function max_img_Callback(hObject, eventdata, handles)


function first_image_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function first_image_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function t_step_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function t_step_Callback(hObject, eventdata, handles)


function run_imEdgeTracker_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% %check if there is a valid image input
% val = get(handles.image_array,'Value');
% if val == 1
%     msgbox('Specify a valid image','Warning','warn');
%     return
% end
% %check if there is a valid working directory input
% val = get(handles.dir_array,'Value');
% if val == 1
%     msgbox('Specify a valid result directory','Warning','warn');
%     return
% end

if handles.parameters_set == 0
    ans = questdlg('You have not specified any parameters! Continue with default settings?','Warning');

    if strcmp(ans,'Yes')
        % get the value from the checkbox
        test_run            = get(handles.test_run,'Value');
        % get the value from the checkbox
        movie               = get(handles.movie,'Value');
        % get the value from the checkbox
        USE_BW_MASK         = get(handles.use_bw_images,'Value');          
        % read the values from the popup menus
        val                 = get(handles.image_array,'Value');
        string_list         = get(handles.image_array,'String');
        firstfilename       = string_list{val}; 
        val                 = get(handles.dir_array,'Value');
        string_list         = get(handles.dir_array,'String');
        hfsmC=findall(0,'Tag','fsmCenter','Name','fsmCenter');
        if ~isempty(hfsmC)
            dir_w = [string_list handles.fsmCenterEdgeDirectory{1} filesep];           
        else
            dir_w = [string_list{val} filesep 'edge' filesep]; 
        end 
        val                 = get(handles.img_depth,'Value');
        string_list         = get(handles.img_depth,'String');
        bit_depth           = string_list{val};  
        if test_run     
            CONTR               = 1;
        else
            CONTR               = 0;           
        end
        PROTRUSION          = get(handles.protrusion,'Value'); 
        FIRST_IMG           = str2num(get(handles.first_image,'String'));
        MAX_IMG             = str2num(get(handles.max_img,'String'));
        T_STEP              = str2num(get(handles.t_step,'String'));
        FILE                = firstfilename;
        RESULTS             = dir_w;
        BIT_DEPTH           = str2num(bit_depth);
        MOVIE               = movie;
    
        %this check should be modified
        if strcmp(FILE,'-')
            msgbox('Specify first a valid image','Warning','warn');
            return
        elseif strcmp(RESULTS,'-')
            msgbox('Specify first a valid directory for the results','Warning','warn');
            return   
        end      

        [img_proccessed, img_edge]=imEdgeTracker(...
            'contr',            CONTR,...
            'file',             FILE,...
            'results',          RESULTS,...
            'protrusion',       PROTRUSION,...
            't_step',           T_STEP,...
            'first_img',        FIRST_IMG,...
            'max_img',          MAX_IMG,...
            'bit_depth',        BIT_DEPTH,... 
            'debug',            handles.parameters.debug,...  
            'movie',            MOVIE,...
            'use_bw_mask',      USE_BW_MASK...
            );
    elseif strcmp(ans,'No')
        return
    elseif strcmp(ans,'Cancel')
        return
    end
else
    %get the value from the checkbox
    test_run            = get(handles.test_run,'Value');
    %get the value from the checkbox
    movie               = get(handles.movie,'Value');
    % get the value from the checkbox
    USE_BW_MASK         = get(handles.use_bw_images,'Value'); 
    %read the values from the popup menus
    val                 = get(handles.image_array,'Value');
    string_list         = get(handles.image_array,'String');
    firstfilename       = string_list{val}; 
    val                 = get(handles.dir_array,'Value');
    string_list         = get(handles.dir_array,'String');
    hfsmC=findall(0,'Tag','fsmCenter','Name','fsmCenter');
    if ~isempty(hfsmC) 
        dir_w = [string_list handles.fsmCenterEdgeDirectory{1} filesep]; 
    else
        dir_w = [string_list{val} filesep 'edge' filesep];
    end
    val                 = get(handles.img_depth,'Value');
    string_list         = get(handles.img_depth,'String');
    bit_depth           = string_list{val}; 
    if test_run
        CONTR               = 1;
    else    
        CONTR               = 0;        
    end
    PROTRUSION          = get(handles.protrusion,'Value');
    MAX_IMG             = str2num(get(handles.max_img,'String'));
    FIRST_IMG           = str2num(get(handles.first_image,'String'));   
    T_STEP              = str2num(get(handles.t_step,'String'));
    FILE                = firstfilename;
    RESULTS             = dir_w;
    BIT_DEPTH           = str2num(bit_depth);   
    MOVIE               = movie;   
    
    if handles.parameters.manual_thresh == 1
        parameters.manual_level = handles.parameters.manual_level;
    else
        parameters.manual_level = -1;
    end
    
    % Test for valid file and directory inputs
    if strcmp(FILE,'-')
        msgbox('Specify first a valid image','Warning','warn');
        return
    elseif strcmp(RESULTS,'-')
        msgbox('Specify first a valid directory for the results','Warning','warn');
        return   
    end 
    
    
    parameters.debug = 0;
    handles.parameters.debug = parameters.debug;
    
    
    [img_proccessed, img_edge]=imEdgeTracker(...
        'contr',            CONTR,...
        'file',             FILE,...
        'results',          RESULTS,...
        't_step',           T_STEP,...
        'first_img',        FIRST_IMG,...       
        'max_img',          MAX_IMG,...
        'pixel',            handles.parameters.pixel,...  
        'time_interval',    handles.parameters.time_interval,...         
        'bit_depth',        BIT_DEPTH,...
        'protrusion',       PROTRUSION,...
        'prot_sampling',    handles.parameters.prot_sampling,... 
        'parenth_l',        handles.parameters.parenth_l,...
        'parenth_r',        handles.parameters.parenth_r,...
        'filter_image',     handles.parameters.filter_image,...
        'img_sigma',        handles.parameters.img_sigma,...
        'erode_dilate',     handles.parameters.erode_dilate,... 
        'f_window',         handles.parameters.f_window,...
        'f_sigma',          handles.parameters.f_sigma,... 
        'tolerance',        handles.parameters.tolerance,...
        'mechanical',       handles.parameters.mechanical,...
        'nearest',          handles.parameters.nearest,...
        'normal',           handles.parameters.normal,...
        'tol',              handles.parameters.tol,...   
        'robust_min',       handles.parameters.robust_min,...  
        'level_set',        handles.parameters.level_set,...          
        'k_S',              handles.parameters.k_s,...
        'k_W',              handles.parameters.k_w,...  
        'cluster',          handles.parameters.cluster,...         
        'cluster_method',   handles.parameters.cluster_method,...        
        'k_cluster',        handles.parameters.k_cluster,...
        'k_min',            handles.parameters.k_min,...   
        'k_max',            handles.parameters.k_max,...
        'cell_mode',        handles.parameters.cell_mode,...
        'debug',            handles.parameters.debug,...    
        'movie',            MOVIE,...
        'use_bw_mask',      USE_BW_MASK,...
        'manual_level',     parameters.manual_level,...
        'orient',           handles.parameters.orient...
        );
end




% --- Executes on button press in protrusion.
function protrusion_Callback(hObject, eventdata, handles)
% guidata(hObject, handles);


function exit_Callback(hObject, eventdata, handles)
delete(handles.figure1);

function help_Callback(hObject, eventdata, handles)
prHelp;


function about_Callback(hObject, eventdata, handles)
prAbout;

function image_array_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function image_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
firstfilename = string_list{val}; 

% change into old image directory if existing
if ~isempty(handles.img_browse_root)
    cd(handles.img_browse_root);
end
if strcmp(firstfilename,'Browse ..')
    [fileName,dirName] = uigetfile('*.tif','Select image');
    if fileName == 0
        return;
    end
    handles.img_browse_root = dirName;
    firstfilename=[dirName,fileName];
    string_list(handles.nr_image_array)   = cellstr(firstfilename);
    string_list(handles.nr_image_array+1) = cellstr('Browse ..');
    
    set(handles.image_array   ,'String', string_list);
    handles.nr_image_array = handles.nr_image_array+1;
end
%cd(getenv('HOME'));

%set project brose root if empty
if isempty(handles.path_browse_root)
    handles.path_browse_root = handles.img_browse_root;
end
guidata(hObject, handles);

function dir_array_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function dir_array_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
dir_w = string_list{val}; 

if strcmp(dir_w,'Browse ..')
    dir_w = uigetdir(handles.path_browse_root,'Select/create directory for results');
    if dir_w == 0
        return;
    end
    handles.path_browse_root = dir_w;
    string_list(handles.nr_dir_array)   = cellstr(dir_w);
    string_list(handles.nr_dir_array+1) = cellstr('Browse ..');
    
    set(handles.dir_array   ,'String', string_list);
    handles.nr_dir_array = handles.nr_dir_array+1;
    guidata(hObject, handles);
end


%check for a parameter file in this directory
fid = fopen([dir_w filesep 'protrusion' filesep 'parameters.dat'],'r');
if fid ~= -1
    ans = questdlg('Parameter file found. Take its values?');
    if strcmp(ans,'Yes')
        parameters = read_parameters([dir_w filesep 'protrusion' filesep 'parameters.dat'], 0);
        %write parameters to tmp disk
        save_parameters([], parameters);
        handles.parameters = parameters;
        handles.parameters_set = 1; 
        
        %set the parameters
        %set bit depth
        bit_depth = parameters.bit_depth;
        if bit_depth == 2^8-1
             set(handles.img_depth,'Value',1);
        elseif bit_depth == 2^10-1
             set(handles.img_depth,'Value',2);
        elseif bit_depth == 2^12-1           
             set(handles.img_depth,'Value',3);
        elseif bit_depth == 2^14-1           
             set(handles.img_depth,'Value',4);
        elseif bit_depth == 2^16-1   
             set(handles.img_depth,'Value',5);            
        end
        %set first image
        set(handles.first_image,'String',num2str(parameters.first_img));
        %set number of images
        set(handles.max_img,'String',num2str(parameters.max_img));
        %set image increment
        set(handles.t_step,'String',num2str(parameters.t_step));
        %set image file
        string_list = get(handles.image_array,'String');
        string_list{end} = parameters.file;
        string_list(handles.nr_image_array+1) = cellstr('Browse ..');
        set(handles.image_array   ,'String', string_list);
        set(handles.image_array   ,'Value', handles.nr_image_array);
        handles.nr_image_array = handles.nr_image_array+1;
        
    end
end
guidata(hObject, handles);



function help_menu_Callback(hObject, eventdata, handles)

function about_menu_Callback(hObject, eventdata, handles)
prAbout;

function helph_menu_Callback(hObject, eventdata, handles)
dir_path_full = which('prPanel');
[pathstr,name,ext,versn] = fileparts(dir_path_full);
web (['file:///', pathstr, filesep, 'prHelpMain.html']); 


function img_depth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function img_depth_Callback(hObject, eventdata, handles)

function test_run_Callback(hObject, eventdata, handles)

function movie_Callback(hObject, eventdata, handles)

function File_Callback(hObject, eventdata, handles)

function use_bw_images_Callback(hObject, eventdata, handles)

function parameters_Callback(hObject, eventdata, handles)



