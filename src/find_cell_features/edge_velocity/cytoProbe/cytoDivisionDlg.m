function varargout = cytoDivisionDlg(varargin)
% CYTODIVISIONDLG M-file for cytoDivisionDlg.fig
%      CYTODIVISIONDLG by itself, creates a new CYTODIVISIONDLG or raises the
%      existing singleton*.
%
%      H = CYTODIVISIONDLG returns the handle to a new CYTODIVISIONDLG or the handle to
%      the existing singleton*.
%
%      CYTODIVISIONDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYTODIVISIONDLG.M with the given input arguments.
%
%      CYTODIVISIONDLG('Property','Value',...) creates a new CYTODIVISIONDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cytoDivisionDlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cytoDivisionDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cytoDivisionDlg

% Last Modified by GUIDE v2.5 24-Aug-2005 16:25:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytoDivisionDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @cytoDivisionDlg_OutputFcn, ...
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

% --- Executes just before cytoDivisionDlg is made visible.
function cytoDivisionDlg_OpeningFcn(hObject, eventdata, handles, varargin)
%set(handles.first_frame,  'Value',varargin{1,1}.first_frame);
%set(handles.roi_bg,      'Value',varargin{1,1}.strech);


handles.output = 'Yes';
guidata(hObject, handles);


if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'title'
          set(hObject, 'Name', varargin{index+1});
         case 'string'
          set(handles.text1, 'String', varargin{index+1});
        end
    end
end

FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end


% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

set(handles.stretch_coeff,'Enable','off');
set(handles.stretch_no,'Value', 1);

% UIWAIT makes cytoDivisionDlg wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = cytoDivisionDlg_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
delete(handles.figure1);

% --- Executes on button press OK
function pushbutton1_Callback(hObject, eventdata, handles)
handles.output.stretch      = get(handles.stretch,'Value');
handles.output.stretch_i    = get(handles.stretch_i,'Value');
handles.output.stretch_c    = get(handles.stretch_c,'Value');
handles.output.stretch_coeff = str2double(get(handles.stretch_coeff,'string'));

handles.output.first_image  = get(handles.first_image,'Value');
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button Chancel
function pushbutton2_Callback(hObject, eventdata, handles)
handles.output.stretch       = get(handles.stretch,'Value');
handles.output.stretch_i     = get(handles.stretch_i,'Value');
handles.output.stretch_c     = get(handles.stretch_c,'Value');
handles.output.stretch_coeff = str2double(get(handles.stretch_coeff,'string'));

handles.output.first_image   = get(handles.first_image,'Value');
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    


function stretch_Callback(hObject, eventdata, handles)
if get(handles.stretch,'Value');
    set(handles.stretch_i,'Value', 0);
    set(handles.stretch_c,'Value', 0);
    set(handles.stretch_no,'Value', 0);
    set(handles.stretch_coeff,'Enable','off');
end

function stretch_i_Callback(hObject, eventdata, handles)
if get(handles.stretch_i,'Value');
    set(handles.stretch,'Value', 0);
    set(handles.stretch_c,'Value', 0);
    set(handles.stretch_no,'Value', 0);
    set(handles.stretch_coeff,'Enable','off');
end

function stretch_c_Callback(hObject, eventdata, handles)
if get(handles.stretch_c,'Value');
    set(handles.stretch,'Value', 0);
    set(handles.stretch_i,'Value', 0);
    set(handles.stretch_no,'Value', 0);
    set(handles.stretch_coeff,'Enable','on');
end


function stretch_no_Callback(hObject, eventdata, handles)
if get(handles.stretch_no,'Value');
    set(handles.stretch,'Value', 0);
    set(handles.stretch_c,'Value', 0);
    set(handles.stretch_i,'Value', 0);
    set(handles.stretch_coeff,'Enable','off');
end



function first_image_Callback(hObject, eventdata, handles)


function stretch_coeff_Callback(hObject, eventdata, handles)
function stretch_coeff_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


