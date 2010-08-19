function varargout = cytoParameterDlg(varargin)
% CYTOPARAMETERDLG M-file for cytoParameterDlg.fig
%      CYTOPARAMETERDLG by itself, creates a new CYTOPARAMETERDLG or raises the
%      existing singleton*.
%
%      H = CYTOPARAMETERDLG returns the handle to a new CYTOPARAMETERDLG or the handle to
%      the existing singleton*.
%
%      CYTOPARAMETERDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYTOPARAMETERDLG.M with the given input arguments.
%
%      CYTOPARAMETERDLG('Property','Value',...) creates a new CYTOPARAMETERDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cytoParameterDlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cytoParameterDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cytoParameterDlg

% Last Modified by GUIDE v2.5 17-Jun-2005 17:16:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytoParameterDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @cytoParameterDlg_OutputFcn, ...
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

% --- Executes just before cytoParameterDlg is made visible.
function cytoParameterDlg_OpeningFcn(hObject, eventdata, handles, varargin)
set(handles.pixel_size,  'String',num2str(varargin{1,1}.pixel_size));
set(handles.pixel_size_c,'Value', varargin{1,1}.pixel_size_c);
set(handles.frame_interval,  'String',num2str(varargin{1,1}.frame_interval));
set(handles.frame_interval_c,'Value', varargin{1,1}.frame_interval_c);

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

% UIWAIT makes cytoParameterDlg wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = cytoParameterDlg_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
delete(handles.figure1);

% --- Executes on button press OK
function pushbutton1_Callback(hObject, eventdata, handles)
handles.output.pixel_size_c         = get(handles.pixel_size_c,'Value');
handles.output.pixel_size           = str2num(get(handles.pixel_size,'String'));
handles.output.frame_interval_c     = get(handles.frame_interval_c,'Value');
handles.output.frame_interval       = str2num(get(handles.frame_interval,'String'));

guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button Chancel
function pushbutton2_Callback(hObject, eventdata, handles)
handles.output.pixel_size_c = get(handles.pixel_size_c,'Value');
handles.output.pixel_size   = str2num(get(handles.pixel_size,'String'));
handles.output.frame_interval_c     = get(handles.frame_interval_c,'Value');
handles.output.frame_interval       = str2num(get(handles.frame_interval,'String'));
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



function pixel_size_Callback(hObject, eventdata, handles)

function pixel_size_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pixel_size_c_Callback(hObject, eventdata, handles)
stat = get(handles.pixel_size_c,'Value');
if stat
    set(handles.pixel_size,'Enable','on');
else
    set(handles.pixel_size,'Enable','off');
end


function frame_interval_Callback(hObject, eventdata, handles)
function frame_interval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function frame_interval_c_Callback(hObject, eventdata, handles)



