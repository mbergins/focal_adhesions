function varargout = cytoRacFRETDlg(varargin)
% CYTORACFRETDLG M-file for cytoRacFRETDlg.fig
%      CYTORACFRETDLG by itself, creates a new CYTORACFRETDLG or raises the
%      existing singleton*.
%
%      H = CYTORACFRETDLG returns the handle to a new CYTORACFRETDLG or the handle to
%      the existing singleton*.
%
%      CYTORACFRETDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYTORACFRETDLG.M with the given input arguments.
%
%      CYTORACFRETDLG('Property','Value',...) creates a new CYTORACFRETDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cytoRacFRETDlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cytoRacFRETDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cytoRacFRETDlg

% Last Modified by GUIDE v2.5 19-Aug-2005 16:16:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytoRacFRETDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @cytoRacFRETDlg_OutputFcn, ...
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

% --- Executes just before cytoRacFRETDlg is made visible.
function cytoRacFRETDlg_OpeningFcn(hObject, eventdata, handles, varargin)
set(handles.CFP_bleedthrough,  'String',num2str(varargin{1,1}));
set(handles.YFP_bleedthrough,  'String',num2str(varargin{1,2}));

% Choose default command line output for cytoRacFRETDlg
handles.output = 'Yes';

% Update handles structure
guidata(hObject, handles);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.
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

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
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
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes cytoRacFRETDlg wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = cytoRacFRETDlg_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in pushbutton1.
function OK_Callback(hObject, eventdata, handles)
background          = get(handles.background, 'Value');
CFP_bleedthrough    = str2double(get(handles.CFP_bleedthrough, 'String'));
YFP_bleedthrough    = str2double(get(handles.YFP_bleedthrough, 'String'));
handles.output.background       = background;
handles.output.CFP_bleedthrough = CFP_bleedthrough;
handles.output.YFP_bleedthrough = YFP_bleedthrough;

guidata(hObject, handles);
uiresume(handles.figure1);


function cancel_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
handles.output = 0;
uiresume(handles.figure1);


function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1);
end


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


function background_Callback(hObject, eventdata, handles)

function CFP_bleedthrough_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function CFP_bleedthrough_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function YFP_bleedthrough_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function YFP_bleedthrough_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


