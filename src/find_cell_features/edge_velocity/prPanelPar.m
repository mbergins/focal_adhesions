function varargout = prPanelPar(varargin)
% PRPANELPAR M-file for prPanelPar.fig
%      PRPANELPAR, by itself, creates a new PRPANELPAR or raises the existing
%      singleton*.
%
%      H = PRPANELPAR returns the handle to a new PRPANELPAR or the handle to
%      the existing singleton*.
%
%      PRPANELPAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRPANELPAR.M with the given input arguments.
%
%      PRPANELPAR('Property','Value',...) creates a new PRPANELPAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prPanelPar_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prPanelPar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prPanelPar

% Last Modified by GUIDE v2.5 24-Aug-2006 10:59:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prPanelPar_OpeningFcn, ...
                   'gui_OutputFcn',  @prPanelPar_OutputFcn, ...
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


% --- Executes just before prPanelPar is made visible.
function prPanelPar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prPanelPar (see VARARGIN)

if nargin > 3
    if strcmp(varargin{1},'parameters')
        parameters = varargin{2};
        %set the fields with the given parameters
        set_values(hObject, eventdata, handles, parameters);
    end
    if strcmp(varargin{3},'status')
        fsmCenterRunning = varargin{4};
    else
        fsmCenterRunning = 0;
    end
end

% fsmCenter is running then take the pixel size and time interval values
% from the fsmCenter panel
if fsmCenterRunning
    set(handles.pixel,          'Enable','off');
    set(handles.time_interval,  'Enable','off');
end

% disable the Level-Set option
set(handles.levelSet,  'Enable','off');

% disable the advanced options
set(handles.levelSet,  'Enable','off');
set(handles.parenth_l,  'Enable','off');
set(handles.parenth_r,  'Enable','off');


% if we use BW images disable the filter options
if parameters.use_bw_images
    set(handles.filter_image, 'Enable','off');
    set(handles.img_sigma,    'Enable','off');
    set(handles.median_f_b, 'Enable','off');
    set(handles.median_f, 'Enable','off');
    set(handles.erode_dilate_button, 'Enable','off');      
    set(handles.erode_dilate, 'Enable','off');   
    % set the manuel threshold to the right values 
    % and disable them
    set(handles.manual_thresh, 'Value',1);
    set(handles.manual_thresh, 'Enable','off');  
    set(handles.manual_level, 'String','0');   
    set(handles.manual_level, 'Enable','off');   
end


guidata(hObject, handles);

set(handles.figure2,'Visible','on');
waitfor(handles.figure2,'Visible');



function prPanelPar_CloseRequestFcn(hObject, eventdata, handles)
% hide the figure here, it will be deleted in the OutputFcn
set(handles.figure2,'Visible','off');

function Ok_Callback(hObject, eventdata, handles)
prPanelPar_CloseRequestFcn(hObject, eventdata, handles);



% --- Outputs from this function are returned to the command line.
function varargout = prPanelPar_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%see http://www.mathworks.com/support/solutions/data/35645.shtml

em              =   get(handles.em,'Value');
if em == 1
    METHOD = 'em';
else
    METHOD = 'kmeans';
end
handles.cluster_method = METHOD;


parameters.pixel           =   str2double(get(handles.pixel,       'String'));
parameters.prot_sampling   =   str2num(get(handles.prot_sampling,  'String'));
parameters.prot_depth      =   0;
parameters.nr_sect         =   0;
parameters.parenth_l       =   str2num(get(handles.parenth_l,      'String'));
parameters.parenth_r       =   str2num(get(handles.parenth_r,      'String'));
parameters.filter_image    =   get(handles.filter_image,           'Value');
parameters.img_sigma       =   str2double(get(handles.img_sigma,   'String'));
parameters.erode_dilate    =   str2num(get(handles.erode_dilate,   'String'));
parameters.median_f        =   get(handles.median_f_b,             'Value');;
parameters.f_window        =   str2num(get(handles.f_window,       'String'));
parameters.f_sigma         =   str2double(get(handles.f_sigma,     'String'));
parameters.tolerance       =   str2double(get(handles.tolerance,   'String'));
parameters.mechanical      =   get(handles.mechanical,             'Value');
parameters.nearest         =   get(handles.nearest,                'Value');
parameters.normal          =   get(handles.normal,                 'Value');
parameters.robust_min      =   get(handles.robust_min,             'Value');
parameters.level_set       =   get(handles.levelSet,               'Value');
parameters.tol             =   str2double(get(handles.tol,         'String'));
parameters.k_s             =   str2double(get(handles.k_s,         'String'));
parameters.k_w             =   str2double(get(handles.k_w,         'String'));
parameters.time_interval   =   str2double(get(handles.time_interval,'String'));
parameters.cluster         =   get(handles.cluster,                'Value');
parameters.cluster_method  =   METHOD;
parameters.k_cluster       =   str2double(get(handles.k_cluster,   'String'));
parameters.k_min           =   str2double(get(handles.k_min,       'String'));
parameters.k_max           =   str2double(get(handles.k_max,       'String'));
parameters.cell_mode       =   str2double(get(handles.cell_modes,  'String'));
parameters.manual_thresh   =   get(handles.manual_thresh,          'Value');
parameters.manual_level    =   str2double(get(handles.manual_level,'String'));

if get(handles.extrema, 'Value');
    if get(handles.orient_left,       'Value');
        parameters.orient  = 0;
    elseif get(handles.orient_right,  'Value');
        parameters.orient  = 1;
    elseif get(handles.orient_up,     'Value');
        parameters.orient  = 2;
    elseif get(handles.orient_low,    'Value');
        parameters.orient  = 3;
    end
elseif get(handles.crosshair, 'Value');  
    if get(handles.orient_left,       'Value');
        parameters.orient  = 4;
    elseif get(handles.orient_right,  'Value');
        parameters.orient  = 5;
    elseif get(handles.orient_up,     'Value');
        parameters.orient  = 6;
    elseif get(handles.orient_low,    'Value');
        parameters.orient  = 7;
    end
elseif get(handles.orient_ellipse,'Value');  
    parameters.orient  = 10;   
end

varargout{1}    =   parameters;
delete(hObject);




function filter_image_Callback(hObject, eventdata, handles)
prot_value = get(hObject,'Value');
if prot_value
    set(handles.img_sigma,'Enable','on');
else
    set(handles.img_sigma,'Enable','off');
end
guidata(hObject, handles);


function pixel_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function pixel_Callback(hObject, eventdata, handles)


function prot_sampling_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function prot_sampling_Callback(hObject, eventdata, handles)


function img_sigma_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function img_sigma_Callback(hObject, eventdata, handles)

function erode_dilate_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function erode_dilate_button_Callback(hObject, eventdata, handles)
prot_value = get(hObject,'Value');
if prot_value
    set(handles.erode_dilate,'Enable','on');
    set(handles.erode_dilate,'String',6);   
else
    set(handles.erode_dilate,'String',0);  
    set(handles.erode_dilate,'Enable','off');
end
guidata(hObject, handles);

function erode_dilate_Callback(hObject, eventdata, handles)

function f_window_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function f_window_Callback(hObject, eventdata, handles)

function f_sigma_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function f_sigma_Callback(hObject, eventdata, handles)

function tolerance_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function tolerance_Callback(hObject, eventdata, handles)


function k_s_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function k_s_Callback(hObject, eventdata, handles)

function k_w_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function k_w_Callback(hObject, eventdata, handles)


function mechanical_Callback(hObject, eventdata, handles)
value = get(hObject,'Value');
if value
    set(handles.k_s,'Enable','on');
    set(handles.k_w,'Enable','on');
    set(handles.robust_min,'Enable','on');    
    set(handles.tol,'Enable','on'); 
    set(handles.nearest,'Value',0);
    set(handles.normal,'Value',0);    
    set(handles.levelSet,'Value',0);    
else
    set(handles.k_s,'Enable','off');
    set(handles.k_w,'Enable','off'); 
end
guidata(hObject, handles);


function nearest_Callback(hObject, eventdata, handles)
value = get(hObject,'Value');
if value
     set(handles.mechanical,'Value',0);
     set(handles.normal,'Value',0);   
     set(handles.levelSet,'Value',0);  
     
     set(handles.robust_min,'Enable','on');    
     set(handles.tol,'Enable','on'); 
end


function normal_Callback(hObject, eventdata, handles)
value = get(hObject,'Value');
if value
     set(handles.mechanical,'Value',0);
     set(handles.nearest,'Value',0);  
     set(handles.levelSet,'Value',0);   
     
     set(handles.robust_min,'Enable','on');    
     set(handles.tol,'Enable','on');      
end


function parenth_l_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function parenth_l_Callback(hObject, eventdata, handles)

function parenth_r_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function parenth_r_Callback(hObject, eventdata, handles)


function median_f_b_Callback(hObject, eventdata, handles)
prot_value = get(hObject,'Value');

if prot_value
    set(handles.median_f,'Enable','on');
    set(handles.median_f,'String',3);   
else
    set(handles.median_f,'String',0);  
    set(handles.median_f,'Enable','off');
end
guidata(hObject, handles);


function median_f_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function median_f_Callback(hObject, eventdata, handles)


function time_interval_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function time_interval_Callback(hObject, eventdata, handles)

function cluster_Callback(hObject, eventdata, handles)
state = get(hObject,'Value');
if state
    set(handles.f_window,'Enable','off');  
    set(handles.f_sigma,'Enable','off');   
    set(handles.manual_thresh,'Value',0);
    set(handles.manual_level,'Enable','off');   
else
    set(handles.f_window,'Enable','on');  
    set(handles.f_sigma,'Enable','on');
    set(handles.manual_level,'Enable','on');       
end



function kmeans_Callback(hObject, eventdata, handles)
state = get(hObject,'Value');
if state
    set(handles.em,'Value',0);  
else
    set(handles.em,'Value',1);  
    set(handles.cell_modes_slider,'Min',1);
    set(handles.cell_modes_slider,'Max',2); 
    set(handles.cell_modes,'Value',get(handles.k_cluster,'Value')); 
end

guidata(hObject, handles);



function em_Callback(hObject, eventdata, handles)
state = get(hObject,'Value');
if state
    set(handles.kmeans,'Value',0);  
else
    set(handles.kmeans,'Value',1);  
    set(handles.cell_modes_slider,'Min',1);
    set(handles.cell_modes_slider,'Max',get(handles.k_cluster,'Value'));
    set(handles.cell_modes,'Value',get(handles.k_max,'Value'));     
end

function k_max_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function k_max_Callback(hObject, eventdata, handles)
set(handles.cell_modes_slider,'Max',str2num(get(handles.k_max,'String')));
set(handles.cell_modes,'Value',str2num(get(handles.k_max,'String')));
guidata(hObject, handles);



function k_cluster_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function k_cluster_Callback(hObject, eventdata, handles)
%set(handles.cell_modes_slider,'Min',1);
set(handles.cell_modes_slider,'Max',str2num(get(handles.k_cluster,'String')));
set(handles.cell_modes,'Value',str2num(get(handles.k_cluster,'String')));
guidata(hObject, handles);



function k_min_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function k_min_Callback(hObject, eventdata, handles)




% --- Executes on slider movement.
function cell_modes_slider_Callback(hObject, eventdata, handles)
% hObject    handle to cell_modes_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.cell_modes,'String',num2str(get(handles.cell_modes_slider,'Value')));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function cell_modes_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Max',3);
set(hObject,'Value',3);
set(hObject,'Min',1);
set(hObject,'SliderStep',[1 1]);
set(hObject,'Value',1);
guidata(hObject, handles);


function cell_modes_Callback(hObject, eventdata, handles)
function cell_modes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function robust_min_Callback(hObject, eventdata, handles)



function tol_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tol_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_Callback(hObject, eventdata, handles)
[fileName,dirName] = uiputfile('*.par','Save as');
if fileName == 0
    return;
end
filePath = [dirName,fileName];

fid = fopen([filePath '.dat'],'w+');
if fid == -1
    error('Could not create parameter file');
end
fclose(fid);


function load_Callback(hObject, eventdata, handles)
[fileName,dirName] = uigetfile('*.par','Open parameter file');
if fileName == 0
    return;
end
filePath = [dirName,fileName];



function set_values(hObject, eventdata, handles, parameters)

set(handles.pixel,          'String',num2str(parameters.pixel));
set(handles.time_interval,  'String',num2str(parameters.time_interval));
set(handles.prot_sampling,  'String',num2str(parameters.prot_sampling));
set(handles.parenth_l,      'String',num2str(parameters.parenth_l));
set(handles.parenth_r,      'String',num2str(parameters.parenth_r));
set(handles.filter_image,   'Value', parameters.filter_image);
set(handles.img_sigma,      'String',num2str(parameters.img_sigma ));
set(handles.erode_dilate,   'String',num2str(parameters.erode_dilate));
set(handles.median_f_b,     'Value', parameters.median_f);
set(handles.f_window,       'String',num2str(parameters.f_window));
set(handles.f_sigma,        'String',num2str(parameters.f_sigma));
set(handles.tolerance,      'String',num2str(parameters.tolerance));
set(handles.mechanical,     'Value', parameters.mechanical);
set(handles.nearest,        'Value', parameters.nearest);
set(handles.normal,         'Value', parameters.normal);
set(handles.robust_min,     'Value', parameters.robust_min);
set(handles.tol,            'String',num2str(parameters.tol));
set(handles.k_s,            'String',num2str(parameters.k_s));
set(handles.k_w,            'String',num2str(parameters.k_w));
set(handles.cluster,        'Value', parameters.cluster);
if strcmp(parameters.cluster_method,'em');
    set(handles.kmeans,'Value',0); 
    set(handles.em,    'Value',1);   
else
    set(handles.kmeans,'Value',1); 
    set(handles.em,    'Value',0);   
end
set(handles.k_cluster,      'String',num2str(parameters.k_cluster));
set(handles.k_min,          'String',num2str(parameters.k_min));
set(handles.k_max,          'String',num2str(parameters.k_max));
if parameters.manual_thresh
    set(handles.manual_thresh,   'Value',1);
    set(handles.manual_level,   'String',num2str(parameters.manual_level));
else
     set(handles.manual_thresh,   'Value',0);   
     set(handles.manual_level,   'String',num2str('0'));   
end
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 0);
set(handles.orient_ellipse,   'Value', 0);
if parameters.orient == 0
    set(handles.extrema,      'Value', 1);    
    set(handles.orient_left,  'Value', 1);
elseif parameters.orient == 1 
    set(handles.extrema,      'Value', 1);       
    set(handles.orient_right, 'Value', 1); 
elseif parameters.orient == 2     
    set(handles.extrema,      'Value', 1);       
    set(handles.orient_up,    'Value', 1);  
elseif parameters.orient == 3    
    set(handles.extrema,      'Value', 1);       
    set(handles.orient_low,   'Value', 1);   
elseif parameters.orient == 4    
    set(handles.crosshair,      'Value', 1);        
    set(handles.orient_left,  'Value', 1);
elseif parameters.orient == 5   
    set(handles.crosshair,      'Value', 1);       
    set(handles.orient_right, 'Value', 1); 
elseif parameters.orient == 6    
    set(handles.crosshair,      'Value', 1);       
    set(handles.orient_up,    'Value', 1);      
elseif parameters.orient == 7    
    set(handles.crosshair,      'Value', 1);       
    set(handles.orient_low,   'Value', 1);   
elseif parameters.orient == 10    
    set(handles.orient_ellipse,  'Value', 1);      
end
set(handles.cell_modes, 'String',num2str(parameters.cell_mode));

%set(handles.cell_modes_slider,'Max',3);
%set(handles.cell_modes_slider,'Min',1);
%set(handles.cell_modes_slider,'SliderStep',[1 1]);

function set_default_values(hObject, eventdata, handles, parameters)
set(handles.pixel,          'String',num2str(67));
set(handles.time_interval,  'String',num2str(10));
set(handles.prot_sampling,  'String',num2str(5));
set(handles.parenth_l,      'String',num2str(0));
set(handles.parenth_r,      'String',num2str(0));
set(handles.filter_image,   'Value', 1);
set(handles.img_sigma,      'String',num2str(0.9));
set(handles.erode_dilate,   'String',num2str(6));
set(handles.median_f_b,     'Value', 1);;
set(handles.f_window,       'String',num2str(5));
set(handles.f_sigma,        'String',num2str(0.1));
set(handles.tolerance,      'String',num2str(20));
set(handles.mechanical,     'Value', 1);
set(handles.nearest,        'Value', 0);
set(handles.normal,         'Value', 0);
set(handles.tol,            'String',num2str(0));
set(handles.robust_min,     'Value', 1);
set(handles.k_s,            'String',num2str(0.01));
set(handles.k_w,            'String',num2str(1));
set(handles.cluster,        'Value', 0);
set(handles.kmeans,         'Value',0); 
set(handles.em,             'Value',1);   
set(handles.k_cluster,      'String',num2str(3));
set(handles.k_min,          'String',num2str(3));
set(handles.k_max,          'String',num2str(3)); 
set(handles.extrema,        'Value', 1);  
set(handles.orient_left,      'Value', 1);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 0);
set(handles.orient_ellipse,   'Value', 0);
set(handles.cell_modes,        'String',num2str(0));


%set(handles.cell_modes_slider,'Value',1);
%set(handles.cell_modes_slider,'Max',3);
%set(handles.cell_modes_slider,'Min',1);
%set(handles.cell_modes_slider,'SliderStep',[1 1]);

function file_Callback(hObject, eventdata, handles)
function exit_Callback(hObject, eventdata, handles)
function context1_Callback(hObject, eventdata, handles)
function open_Callback(hObject, eventdata, handles)
function Untitled_2_Callback(hObject, eventdata, handles)


function manual_level_Callback(hObject, eventdata, handles)

function manual_level_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function manual_thresh_Callback(hObject, eventdata, handles)
if get(handles.manual_thresh,'Value')
    set(handles.cluster,'Value',0);
    set(handles.manual_level,'Enable','on');    
end


function levelSet_Callback(hObject, eventdata, handles)

if get(handles.levelSet,'Value')
    set(handles.mechanical,'Value',0);
    set(handles.nearest,'Value',0);  
    set(handles.normal,'Value',0); 
    set(handles.robust_min,'Enable','off');      
    set(handles.tol,'Enable','off'); 
else
    set(handles.mechanical,'Value',1);
    set(handles.nearest,'Value',0);  
    set(handles.normal,'Value',0);     
    set(handles.robust_min,'Enable','on');    
    set(handles.tol,'Enable','on');   
end

function orient_ellipse_Callback(hObject, eventdata, handles)
set(handles.orient_left,  'Enable', 'off');
set(handles.orient_right, 'Enable', 'off');
set(handles.orient_up,    'Enable', 'off');
set(handles.orient_low,   'Enable', 'off');

set(handles.orient_ellipse,'Value', 1);
set(handles.crosshair,     'Value', 0);
set(handles.extrema,     'Value', 0);

function extrema_Callback(hObject, eventdata, handles)
set(handles.orient_left,  'Enable', 'on');
set(handles.orient_right, 'Enable', 'on');
set(handles.orient_up,    'Enable', 'on');
set(handles.orient_low,   'Enable', 'on');

set(handles.orient_ellipse, 'Value', 0);
set(handles.crosshair,      'Value', 0);

function crosshair_Callback(hObject, eventdata, handles)
set(handles.orient_left,  'Enable', 'on');
set(handles.orient_right, 'Enable', 'on');
set(handles.orient_up,    'Enable', 'on');
set(handles.orient_low,   'Enable', 'on');

set(handles.orient_ellipse,'Value', 0);
set(handles.extrema,       'Value', 0);


function orient_left_Callback(hObject, eventdata, handles)
set(handles.orient_left,      'Value', 1);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 0);
set(handles.orient_ellipse,   'Value', 0);

function orient_right_Callback(hObject, eventdata, handles)
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 1);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 0);
set(handles.orient_ellipse,   'Value', 0);

function orient_up_Callback(hObject, eventdata, handles)
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 1);
set(handles.orient_low,       'Value', 0);
set(handles.orient_ellipse,   'Value', 0);

function orient_low_Callback(hObject, eventdata, handles)
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 1);
set(handles.orient_ellipse,   'Value', 0);

function orient_ellipse_Callback(hObject, eventdata, handles)
set(handles.orient_left,      'Value', 0);
set(handles.orient_right,     'Value', 0);
set(handles.orient_up,        'Value', 0);
set(handles.orient_low,       'Value', 0);
set(handles.orient_ellipse,   'Value', 1);




