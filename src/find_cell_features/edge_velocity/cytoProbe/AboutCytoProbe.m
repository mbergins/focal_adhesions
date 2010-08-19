function varargout = AboutCytoProbe(varargin)
% ABOUTCYTOPROBE M-file for AboutCytoProbe.fig
%      ABOUTCYTOPROBE, by itself, creates a new ABOUTCYTOPROBE or raises the existing
%      singleton*.
%
%      H = ABOUTCYTOPROBE returns the handle to a new ABOUTCYTOPROBE or the handle to
%      the existing singleton*.
%
%      ABOUTCYTOPROBE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABOUTCYTOPROBE.M with the given input arguments.
%
%      ABOUTCYTOPROBE('Property','Value',...) creates a new ABOUTCYTOPROBE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AboutCytoProbe_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AboutCytoProbe_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AboutCytoProbe

% Last Modified by GUIDE v2.5 09-Jan-2004 14:07:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AboutCytoProbe_OpeningFcn, ...
                   'gui_OutputFcn',  @AboutCytoProbe_OutputFcn, ...
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


% --- Executes just before AboutCytoProbe is made visible.
function AboutCytoProbe_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AboutCytoProbe (see VARARGIN)


% Choose default command line output for AboutCytoProbe
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AboutCytoProbe wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AboutCytoProbe_OutputFcn(hObject, eventdata, handles)
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


