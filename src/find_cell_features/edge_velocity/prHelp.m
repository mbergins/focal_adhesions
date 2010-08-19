function varargout = prHelp(varargin)
% PRHELP M-file for prHelp.fig
%      PRHELP, by itself, creates a new PRHELP or raises the existing
%      singleton*.
%
%      H = PRHELP returns the handle to a new PRHELP or the handle to
%      the existing singleton*.
%
%      PRHELP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRHELP.M with the given input arguments.
%
%      PRHELP('Property','Value',...) creates a new PRHELP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prHelp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prHelp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prHelp

% Last Modified by GUIDE v2.5 20-Jan-2004 14:42:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prHelp_OpeningFcn, ...
                   'gui_OutputFcn',  @prHelp_OutputFcn, ...
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


% --- Executes just before prHelp is made visible.
function prHelp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prHelp (see VARARGIN)

handles.ind = 1;

string = {'The program (the main function is imEdgeTracker) extracts the boundary', 
'of an object and can determine the advancement of this boundary (the',
'so-called protrusion) in time series images.',
' ',
'The following inrformation can be calculated with the code:'
'- object edge [pixel coordinates],'
'- object image mask [binary image],'
'- object edge band [binary image],'
'- object edge band segments [pixel coordinates],'
'- edge normals [vectors],'
'- edge protrusion [vectors],'
' ',
'The code works as following:',
'The image i_s segmented based on pixel intensity values. The',
'threshold value for this is automaticly determined based on', 
'an analysis of the image histogram.', 
'The segmented (binary) image is median filtered and cleaned up.', 
'Attention, if the cell covers only a small part of the image',
'unexpected results can occure since a size criterion is used',
'to identify the cell body.',
'The cell edge is defined as a 8-nghb. connected pixel object. For',
'the further analysis the edge is approximated with a smoothing',
'spline. The edge protrusion (displacement) is calculated using a',
'mechanical model. This requires a initial solution calculated by',
'the nearest model. A normal model is also implemented but not',
'used in Version 1.0',
'After the stack processing some data analyis is performed.', 
'As a control output the variable "img_proccessed" is given.',
'The values are: 0 if the image is not processed (because',
't_step is not 1), 1 for a successful process and -1 when', 
'the atempt failed.'};    
      
[outstring,newpos] = textwrap(handles.help_text,string);
set(handles.help_text,'String',outstring)


% Choose default command line output for prHelp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prHelp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prHelp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
delete(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)

if handles.ind == 1
    handles.ind = 2;
    string = {'The program (the main function is imEdgeTracker) extracts the boundary', 
        'of an object and can determine the advancement of this boundary (the',
        'so-called protrusion) in time series images.',
        ' ',
        'The following inrformation can be calculated with the code:'
        '- object edge [pixel coordinates],'
        '- object image mask [binary image],'
        '- object edge band [binary image],'
        '- object edge band segments [pixel coordinates],'
        '- edge normals [vectors],'
        '- edge protrusion [vectors],'
        ' ',
        'The code works as following:',
        'The image is segmented based on pixel intensity values. The',
        'threshold value for this is automaticly determined based on', 
        'an analysis of the image histogram.', 
        'The segmented (binary) image is median filtered and cleaned up.', 
        'Attention, if the cell covers only a small part of the image',
        'unexpected results can occure since a size criterion is used',
        'to identify the cell body.',
        'The cell edge is defined as a 8-nghb. connected pixel object. For',
        'the further analysis the edge is approximated with a smoothing',
        'spline. The edge protrusion (displacement) is calculated using a',
        'mechanical model. This requires a initial solution calculated by',
        'the nearest model. A normal model is also implemented but not',
        'used in Version 1.0',
        'After the stack processing some data analyis is performed.', 
        'As a control output the variable "img_proccessed" is given.',
        'The values are: 0 if the image is not processed (because',
        't_step is not 1), 1 for a successful process and -1 when', 
        'the atempt failed.'};    
elseif handles.ind == 2;
    string = {'The data structure of the normals is:',
         '(this 6 lines are repeated for each image)',
         'number of normals',
         'spline parameters of the normal vector position ',        
         'normal vector x-position',
         'normal vector y-position',  
         'normal vector x-dir',
         'normal vector y-dir'
         '',
         '',        
         'The data structure of the protrusion vectors is:',
         '(this 6 lines are repeated for each image)',         
         'number of protrusion vectors',
         'spline parameters of the protrusion vectors position ',        
         'protrusion vectors x-position',
         'protrusion vectors y-position',  
         'protrusion vectors x-dir',
         'protrusion vectors y-dir' 
         '',
         '',  
         'the spline that is calculated for the object boundary at each time step',
         'has parameters running from [1..number of edge pixels]',
     };
     
     
     
    handles.ind = 3;
elseif handles.ind == 3;
    string = {'Just wait more'};   
    handles.ind = 1;  
end

[outstring,newpos] = textwrap(handles.help_text,string);
set(handles.help_text,'String',outstring)
guidata(hObject, handles);

% --- Executes on button press in perv.
function perv_Callback(hObject, eventdata, handles)
if handles.ind == 1
    handles.ind = 2;
    next_Callback(hObject, eventdata, handles);
elseif handles.ind == 2
    handles.ind = 3;
    next_Callback(hObject, eventdata, handles);
elseif handles.ind == 3
    handles.ind = 1;
    next_Callback(hObject, eventdata, handles);
end 