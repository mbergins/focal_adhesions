function varargout = cytoProbe(varargin)
% CYTOPROBE M-file for cytoProbe.fig
%      CYTOPROBE, by itself, creates a new CYTOPROBE or raises the existing
%      singleton*.
%
%      H = CYTOPROBE returns the handle to a new CYTOPROBE or the handle to
%      the existing singleton*.
%
%      CYTOPROBE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYTOPROBE.M with the given input arguments.
%
%      CYTOPROBE('Property','Value',...) creates a new CYTOPROBE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cytoProbe_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cytoProbe_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help cytoProbe

% Last Modified by GUIDE v2.5 03-Aug-2006 14:09:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytoProbe_OpeningFcn, ...
                   'gui_OutputFcn',  @cytoProbe_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cytoProbe is made visible.
function cytoProbe_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cytoProbe (see VARARGIN)

% Choose default command line output for cytoProbe
handles.output = hObject;

% UIWAIT makes cytoProbe wait for user response (see UIRESUME)
% uiwait(handles.cytoProbeFigure);

% fill the image 1 drop down menu
image1_array = cell(1,1);
image1_array(1) = cellstr('-- Choose an image --');
image1_array(2) = cellstr('Browse ..');
set(handles.image1   ,'String', image1_array);
handles.nr_image1_array = 2;

% fill the image 2 drop down menu
image2_array = cell(1,1);
image2_array(1) = cellstr('-- Choose an image --');
image2_array(2) = cellstr('Browse ..');
set(handles.image2   ,'String', image2_array);
handles.nr_image2_array = 2;

% fill the vector drop down menu
vector_array = cell(1,1);
vector_array(1) = cellstr('-- Choose an scaling vector --');
vector_array(2) = cellstr('Browse ..');
set(handles.sVector   ,'String', vector_array);
handles.nr_vector_array = 2;


set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');

% fill the initial path 
handles.path_browse_root = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some specific handles
%FRET tool
handles.CFP_bleedthrough = 0.395;
handles.YFP_bleedthrough = 0.06;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Update handles structure
guidata(hObject, handles);
            
% --- Outputs from this function are returned to the command line.
function varargout = cytoProbe_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function image1_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
dir_w = string_list{val}; 
if strcmp(dir_w,'Browse ..')
    [image1_name, image1_path_name] = uigetfile({'*.tif';'*.tiff';'*.TIF';'*.bmp';'*.mat'},...
                                    handles.path_browse_root,'Select image 1');
    if image1_name == 0
        return;
    end
    handles.path_browse_root = image1_path_name;
    string_list(handles.nr_image1_array)   = cellstr([image1_path_name image1_name]);
    string_list(handles.nr_image1_array+1) = cellstr('Browse ..');
    
    set(handles.image1   ,'String', string_list);
    handles.nr_image1_array = handles.nr_image1_array+1;
    guidata(hObject, handles);
end
cd(handles.path_browse_root);

function image1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function image2_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
dir_w = string_list{val}; 
if strcmp(dir_w,'Browse ..')
    [image2_name, image2_path_name] = uigetfile({'*.tif';'*.tiff';'*.TIF';'*.bmp';'*.mat'},...
                                       handles.path_browse_root,'Select image 2');
    if image2_name == 0
        return;
    end
    handles.path_browse_root = image2_path_name;
    string_list(handles.nr_image2_array)   = cellstr([image2_path_name image2_name]);
    string_list(handles.nr_image2_array+1) = cellstr('Browse ..');
    
    set(handles.image2   ,'String', string_list);
    handles.nr_image2_array = handles.nr_image2_array+1;
    guidata(hObject, handles);
end


function image2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sVector_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
string_list = get(hObject,'String');
dir_w = string_list{val}; 
if strcmp(dir_w,'Browse ..')
    [vector_name, vector_path_name] = uigetfile('*.mat', handles.path_browse_root,'Scaling vector');
    if vector_name == 0
        return;
    end
    handles.path_browse_root = vector_path_name;
    string_list(handles.nr_vector_array)   = cellstr([vector_path_name vector_name]);
    string_list(handles.nr_vector_array+1) = cellstr('Browse ..');
    
    set(handles.sVector   ,'String', string_list);
    handles.nr_vector_array = handles.nr_vector_array+1;
    guidata(hObject, handles);
end



function sVector_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function run_Callback(hObject, eventdata, handles)

% figure out what type of job to do
racFret         = get(handles.rac(1),'Check');
multiply        = get(handles.multiply(1),'Check');
divide          = get(handles.divide(1),'Check');
bleach          = get(handles.get_bleach(1),'Check');
sMultiply       = get(handles.img_vec_mult(1),'Check');
correlate       = get(handles.correlate(1),'Check');
activity_map    = get(handles.activity_from_edge(1),'Check');
inverse         = get(handles.inverse(1),'Check');
subtract        = get(handles.subtract(1),'Check');
subtract_v      = get(handles.subtract_v(1),'Check');
segmentation_a  = get(handles.automatic(1),'Check');
segmentation_m  = get(handles.manual(1),'Check');
bleach_correct  = get(handles.bleach_correct(1),'Check');
curvature       = get(handles.curvature(1),'Check');
subtract_bg     = get(handles.subtract_bg(1),'Check');
translate       = get(handles.translate(1),'Check');
get_translation = get(handles.get_translation(1),'Check');
make_movie      = get(handles.make_movie(1),'Check');
average_image   = get(handles.average_image(1),'Check');
overlay_pixels  = get(handles.overlay(1),'Check');
average_vector_field  = get(handles.average_vector_field(1),'Check');
track_analysis  = get(handles.track_analysis(1),'Check');
mpm_to_vec      = get(handles.mpm_to_vec(1),'Check');
scores_to_vec   = get(handles.scores_to_vec(1),'Check');
analyse_prot    = get(handles.analyse_prot(1),'Check');
pka    = get(handles.pka(1),'Check');


% get how to read the data stacks
read_mode_str      = get(handles.read_mode_s(1),'Check');
if strcmp(read_mode_str,'on')
    read_mode = 0;
else
    read_mode = 1;
end

%%%%%%%%% check for the paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get first image path
val                 = get(handles.image1,'Value');
string_list         = get(handles.image1,'String');
image1_path         = string_list{val};

if length(string_list) == 2
    image1_path_exist = 0;
elseif length(string_list) > 2 & (val==length(string_list) | val==1)
    image1_path_exist = 0;
else
    image1_path_exist = 1;   
end

% get second image path
val                 = get(handles.image2,'Value');
string_list         = get(handles.image2,'String');
image2_path         = string_list{val};

if length(string_list) == 2
    image2_path_exist = 0;    
elseif length(string_list) > 2 & (val==length(string_list) | val==1)
    image2_path_exist = 0;     
else
    image2_path_exist = 1;   
end
    
% get vector path
val                 = get(handles.sVector,'Value');
string_list         = get(handles.sVector,'String');
vector_path         = string_list{val};

if length(string_list) == 2
    vector_path_exist = 0;    
elseif length(string_list) > 2 & (val==length(string_list) | val==1)
    vector_path_exist = 0;     
else
    vector_path_exist = 1;
end

% get factor
factor = str2double(get(handles.factor,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the parameter 
if isfield(handles, 'cytoParameters')
    data.parameters = handles.cytoParameters;
else
    data.parameters.pixel_size_c = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the result path  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if image1_path_exist
    [result_path] = fileparts(image1_path);
    while ~strcmp(result_path(end), filesep) 
        result_path(end) = [];
    end
else
    return
end

if  strcmp(multiply,'on')
    if image1_path_exist & image2_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';
        data.result = result_path;     
        data.factor = factor;        
        data.read_mode = read_mode;
        imOperateImgStack(data, 1);
    else
        h = warndlg('Specify the images','Warning');
    end

elseif strcmp(divide ,'on')   
    if image1_path_exist & image2_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;  
        data.factor = factor;        
        data.read_mode = read_mode;
        imOperateImgStack(data, 2);
    else
        h = warndlg('Specify the images','Warning');    
    end
    
elseif strcmp(sMultiply,'on')
    if image1_path_exist & vector_path_exist
        data.image1_path = image1_path;
        data.image2_path = '';        
        data.sVector_path = vector_path;   
        data.result = result_path; 
        data.factor = factor;        
        imOperateImgStack(data, 3);
    else
        h = warndlg('Specify the images and vectors','Warning');    
    end

elseif strcmp(bleach,'on')
    if image1_path_exist
        data.image1_path = image1_path;
        data.image2_path = '';        
        data.sVector_path = '';     
        data.result = result_path;    
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,4);
    else
        h = warndlg('Specify the images','Warning');    
    end
    
elseif strcmp(racFret,'on')
    if image1_path_exist  & image2_path_exist & factor
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.image3_path = vector_path;     
        data.result = result_path;       
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,5);
    else
        h = warndlg('Specify the images','Warning');    
    end

elseif strcmp(activity_map,'on')
    if image1_path_exist  & image2_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,6);
    else
        h = warndlg('Specify the images','Warning');    
    end
elseif strcmp(correlate,'on')
    if image1_path_exist  & image2_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,7);
    else
        h = warndlg('Specify the images','Warning');    
    end
elseif strcmp(inverse,'on')
    if image1_path_exist
        data.image1_path = image1_path;
        data.image2_path = '';        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,8);
    else
        h = warndlg('Specify the images','Warning');    
    end     
elseif strcmp(subtract,'on')
    if image1_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,9);
    else
        h = warndlg('Specify the images','Warning');    
    end    
elseif strcmp(subtract_v,'on')
    if image1_path_exist & image2_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,10);
    else
        h = warndlg('Specify the images','Warning');    
    end 
elseif strcmp(segmentation_a,'on')
    if image1_path_exist
        data.image1_path = image1_path;
        data.image2_path = '';        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,11);       
    else
        h = warndlg('Specify the images','Warning');    
    end     
elseif strcmp(segmentation_m,'on')
    if image1_path_exist & image2_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = factor;
        data.read_mode = read_mode;
        imOperateImgStack(data,12);
    else
        h = warndlg('Specify the images','Warning');    
    end   
elseif strcmp(bleach_correct,'on')
    if image1_path_exist & image2_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,13);
    else
        h = warndlg('Specify the images','Warning');    
    end       
elseif strcmp(curvature,'on')
    if image1_path_exist
        data.image1_path = image1_path;
        data.image2_path = '';        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,14);
    else
        h = warndlg('Specify the images','Warning');    
    end           
elseif strcmp(subtract_bg,'on')
    if image1_path_exist
        data.image1_path = image1_path;
        data.image2_path = image2_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,15);
    else
        h = warndlg('Specify the images','Warning');    
    end        
elseif strcmp(translate,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,16);
    else
        h = warndlg('Specify the images','Warning');    
    end     
elseif strcmp(get_translation,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = image2_path;
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,17);
    else
        h = warndlg('Specify the images','Warning');    
    end         
elseif strcmp(make_movie,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = '';
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,18);
    else
        h = warndlg('Specify the images','Warning');    
    end      
elseif strcmp(average_image,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = '';
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,19);
    else
        h = warndlg('Specify the images','Warning');    
    end     
elseif strcmp(overlay_pixels,'on')
    if image1_path_exist & image2_path_exist
        data.image1_path = image1_path;        
        data.image2_path = image2_path;
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,20);
    else
        h = warndlg('Specify the images and pixels list','Warning');    
    end       
elseif strcmp(average_vector_field,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = '';
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,21);        
    else
        h = warndlg('Specify the images and pixels list','Warning');    
    end    
elseif strcmp(mpm_to_vec,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = '';
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,22);        
    else
        h = warndlg('Specify the MPM matrix','Warning');    
    end   
elseif strcmp(track_analysis,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = '';
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,23);        
    else
        h = warndlg('Specify MPM matrix','Warning');    
    end  
elseif strcmp(scores_to_vec,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = '';
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,24);        
    else
        h = warndlg('Specify SCORES file','Warning');    
    end       
elseif strcmp(analyse_prot,'on')
    if image1_path_exist
        data.image1_path = image1_path;        
        data.image2_path = '';
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,25);        
    else
        h = warndlg('Specify protrusion file','Warning');    
    end      
elseif strcmp(pka,'on')
    if image1_path_exist
        data.image1_path =  image1_path;        
        data.image2_path =  image2_path;
        data.sVector_path = '';     
        data.result = result_path;       
        data.factor = '';
        data.read_mode = read_mode;
        imOperateImgStack(data,26);        
    else
        h = warndlg('Specify protrusion file','Warning');    
    end        
else
    return;
end


function file_Callback(hObject, eventdata, handles)
function Untitled_1_Callback(hObject, eventdata, handles)
function help_Callback(hObject, eventdata, handles)
function about_Callback(hObject, eventdata, handles)
AboutCytoProbe;

function helphelp_Callback(hObject, eventdata, handles)
dir_path_full = which('cytoProbe');
[pathstr,name,ext,versn] = fileparts(dir_path_full);
web (['file:///', pathstr, filesep, 'HelpCytoProbe.html']); 



function factor_Callback(hObject, eventdata, handles)
function factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
function vector_field_operations_Callback(hObject, eventdata, handles)
function correlate_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 1);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Vector field 1');
set(handles.static_2,'String','Vector field 2');

set(handles.info,'String','Correlates Vector field 1 with vector field 2');

function subtract_v_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 2);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Vector field 1');
set(handles.static_2,'String','Vector field 2');

set(handles.info,'String',['Subtracts Vector field 2 from vector field 1. '...
    '   First output is the difference in speed. Second output is the difference vector']);
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
function image_operations_Callback(hObject, eventdata, handles)
function inverse_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 3);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image 1');
set(handles.static_2,'String','Image 2');

set(handles.info,'String','Inverse of Image 1');

function subtract_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 4);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image 1');
set(handles.static_2,'String','Image 2');

set(handles.info,'String','Subtract Image 2 from Image 1. Values < 0 are set to zero.');


function multiply_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 5);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image 1');
set(handles.static_2,'String','Image 2');

set(handles.info,'String','Multiplies image 1 by Image 2');

% --------------------------------------------------------------------
function divide_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 6);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');

set(handles.static_1,'String','Image 1');
set(handles.static_2,'String','Image 2');
set(handles.info,'String','Divides image 1 by Image 2.');

function get_bleach_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 7);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');

set(handles.static_1,'String','Image 1');
set(handles.static_2,'String','Image 2');
set(handles.info,'String','Gets bleach coefficient from image 1 by quadratic approxiamation of the mean pixel intensity');

function img_vec_mult_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 8);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','on');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image 1');
set(handles.static_2,'String','Image 2');
set(handles.info,'String','Multiplies Image 1 element-wise by Vector');
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
function fret_Callback(hObject, eventdata, handles)
function rac_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 9);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','on');
set(handles.factor,'Enable','off');

set(handles.static_1,'String','FRET image');
set(handles.static_2,'String','CFP image');
set(handles.vector,'String','YFP image');
set(handles.info,'String','If there is no bleedthrough from the YFP channel, no image has to be specified');


%**************************************************************************
%**************************************************************************
% --------------------------------------------------------------------
function special_Callback(hObject, eventdata, handles)
function activity_from_edge_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 10);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');

set(handles.static_1,'String','Image, vector field, scrores');
set(handles.static_2,'String','BW cell mask');
set(handles.info,'String','Gets variable as a function from edge distance. Variable 1 can be an image or a vector field. Image 2: BW mask');

%**************************************************************************
%**************************************************************************



function settings_Callback(hObject, eventdata, handles)

function read_mode_s_Callback(hObject, eventdata, handles)
if strcmp(get(handles.read_mode_s,'Check'),'off')
    set(handles.read_mode_s(1),'Check','on')
else
    set(handles.read_mode_s(1),'Check','off') 
end


% --------------------------------------------------------------------
function parameters_Callback(hObject, eventdata, handles)
if isfield(handles, 'cytoParameters')
    handles.cytoParameters = cytoParameterDlg(handles.cytoParameters);
else
    cytoParameters.pixel_size_c         = 0;
    cytoParameters.pixel_size           = 1;
    cytoParameters.frame_interval_c     = 0;
    cytoParameters.frame_interval       = 1;    
    handles.cytoParameters = cytoParameterDlg(cytoParameters);
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function segmentation_Callback(hObject, eventdata, handles)



function automatic_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 11);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image');
set(handles.static_2,'String','');

set(handles.info,'String','Image segmentation by automatic threshold estimation');

function manual_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 12);

set(handles.info,'String','Image segmentation by manual threshold. Does not work yet');


% --------------------------------------------------------------------
function bleach_correct_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 13);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image');
set(handles.static_2,'String','BW Mask');

set(handles.info,'String','Bleach correction with exponential approx. of the cell intensity');


% --------------------------------------------------------------------
function curvature_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 14);

set(handles.info,'String','Get curvature from splines');


% --------------------------------------------------------------------
function subtract_bg_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 15);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image');
set(handles.static_2,'String','BW Mask');

set(handles.info,'String','Subtract av value from image obtained from region given by mask or by user defined ROI');


% --------------------------------------------------------------------
function translate_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 16);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image');
set(handles.static_2,'String',' ');

set(handles.info,'String','Translates the Image in the x and y direction. (x: Pointing right, y: Pointing down)');


% --------------------------------------------------------------------
function get_translation_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 17);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Ref image');
set(handles.static_2,'String','Trans image');

set(handles.info,'String','Finds the translation and rotation between image 1 and image 2 based on a cross-correlation alg. (Code by Ge Yang)');


% --------------------------------------------------------------------
function make_movie_Callback(hObject, eventdata, handles)

cyto_set_check(handles, 18);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image sequence (tif)');
set(handles.static_2,'String','');

set(handles.info,'String','Generates avi movie from tif image stack with colormap of choice');


% --------------------------------------------------------------------
function average_image_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 19);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image sequence (tif)');
set(handles.static_2,'String','');

set(handles.info,'String','Calculates average image from image stack');


% --------------------------------------------------------------------
function overlay_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 20);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Image');
set(handles.static_2,'String','');

set(handles.info,'String','Overlays pixels on image');


% --------------------------------------------------------------------
function average_vector_field_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 21);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Vector field');
set(handles.static_2,'String','');

set(handles.info,'String','Averages vector field in ROI');


function Untitled_2_Callback(hObject, eventdata, handles)


function mpm_to_vec_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 22);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','MPM matrix');
set(handles.static_2,'String','');

set(handles.info,'String','Converts MPM matrix into a file stack. Each file represents a time step');


function track_analysis_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 23);

set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','MPM matrix');
set(handles.static_2,'String','');

set(handles.info,'String','Analyses speckle track speed and life time given a MPM matrix.');


function scores_to_vec_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 24);
set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','SCORES matrix');
set(handles.static_2,'String','');

set(handles.info,'String','Converts SCORES matrix into a stack of files with scores for each time step.');


function analyse_prot_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 25);
set(handles.image1,'Enable','on');
set(handles.image2,'Enable','off');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','Protrusion matrix');
set(handles.static_2,'String','');

set(handles.info,'String','Analyses protrusion matrix obtained be the edge tracking tool');


function pka_Callback(hObject, eventdata, handles)
cyto_set_check(handles, 26);
set(handles.image1,'Enable','on');
set(handles.image2,'Enable','on');
set(handles.sVector,'Enable','off');
set(handles.factor,'Enable','off');
set(handles.static_1,'String','CFP-YFP image');
set(handles.static_2,'String','Flatfield image');
%set(handles.vector,'String','Flatfield image');

set(handles.info,'String','Generates ratio image from YFP CFP channels on the same image.');

