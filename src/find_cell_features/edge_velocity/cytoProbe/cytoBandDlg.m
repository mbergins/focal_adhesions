function varargout = cytoBandDlg(varargin)
% CYTOBANDDLG M-file for cytoBandDlg.fig
%      CYTOBANDDLG by itself, creates a new CYTOBANDDLG or raises the
%      existing singleton*.
%
%      H = CYTOBANDDLG returns the handle to a new CYTOBANDDLG or the handle to
%      the existing singleton*.
%
%      CYTOBANDDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYTOBANDDLG.M with the given input arguments.
%
%      CYTOBANDDLG('Property','Value',...) creates a new CYTOBANDDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cytoBandDlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cytoBandDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cytoBandDlg

% Last Modified by GUIDE v2.5 27-Mar-2006 11:35:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytoBandDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @cytoBandDlg_OutputFcn, ...
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

% --- Executes just before cytoBandDlg is made visible.
function cytoBandDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cytoBandDlg (see VARARGIN)

% Choose default command line output for cytoBandDlg
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

% UIWAIT makes cytoBandDlg wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = cytoBandDlg_OutputFcn(hObject, eventdata, handles)
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
band_width   = str2num(get(handles.band_width, 'String'));
background   = get(handles.background, 'Value');
roi          = get(handles.roi, 'Value');
scores       = get(handles.scores, 'Value');
vectors      = get(handles.vectors, 'Value');
images       = get(handles.images, 'Value');
nr_user_bands= str2num(get(handles.nr_user_bands, 'String'));
normalise_beg= get(handles.normalise_beg, 'Value');
normalise_end= get(handles.normalise_end, 'Value');
if normalise_beg   
    handles.output.normalise = 1;
elseif normalise_end
    handles.output.normalise = 2;
else
    handles.output.normalise = 0;  
end

handles.output.band_width = band_width;
handles.output.background = background;
handles.output.scores = scores;
handles.output.vectors = vectors;
handles.output.images = images;
handles.output.roi = roi;
handles.output.nr_user_bands = nr_user_bands;

pos_and_neg = get(handles.pos_and_neg, 'Value');
pos = get(handles.pos, 'Value');
neg = get(handles.neg, 'Value');
if pos_and_neg   
    handles.output.scores_values = 0;
elseif pos
    handles.output.scores_values = 1;
else
    handles.output.scores_values = 2; 
end

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


function band_width_Callback(hObject, eventdata, handles)


function band_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function background_Callback(hObject, eventdata, handles)


function nr_user_bands_Callback(hObject, eventdata, handles)
function nr_user_bands_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function normalise_beg_Callback(hObject, eventdata, handles)
if get(handles.normalise_beg, 'Value');
   set(handles.normalise_end, 'Value',0);
end

function normalise_end_Callback(hObject, eventdata, handles)
if get(handles.normalise_end, 'Value');
    set(handles.normalise_beg, 'Value',0);
end


function pos_and_neg_Callback(hObject, eventdata, handles)
if get(handles.pos_and_neg, 'Value');
    set(handles.pos, 'Value',0);
    set(handles.neg, 'Value',0);
end

function pos_Callback(hObject, eventdata, handles)
if get(handles.pos, 'Value');
    set(handles.pos_and_neg, 'Value',0);
    set(handles.neg, 'Value',0);
end

function neg_Callback(hObject, eventdata, handles)
if get(handles.neg, 'Value');
    set(handles.pos, 'Value',0);
    set(handles.pos_and_neg, 'Value',0);
end

function roi_Callback(hObject, eventdata, handles)


function scores_Callback(hObject, eventdata, handles)
if get(handles.scores, 'Value');
    set(handles.neg, 'Enable','on');
    set(handles.pos, 'Enable','on');
    set(handles.pos_and_neg, 'Enable','on');
    
    set(handles.images, 'Value',0);
    set(handles.vectors, 'Value',0);
    set(handles.background, 'Enable','off');
else
    set(handles.neg, 'Enable','off');
    set(handles.pos, 'Enable','off');
    set(handles.pos_and_neg, 'Enable','off');
end

function images_Callback(hObject, eventdata, handles)
if get(handles.images, 'Value');
    set(handles.scores,'Value',0);
    set(handles.vectors, 'Value',0);
    
    set(handles.background, 'Enable','on');
    
    set(handles.neg, 'Enable','off');
    set(handles.pos, 'Enable','off');
    set(handles.pos_and_neg, 'Enable','off');    
end


function vectors_Callback(hObject, eventdata, handles)
if get(handles.vectors, 'Value');
    set(handles.scores,'Value',0);
    set(handles.images,'Value',0);
    
    set(handles.neg, 'Enable','off');
    set(handles.pos, 'Enable','off');
    set(handles.pos_and_neg, 'Enable','off');
    
    set(handles.background, 'Enable','off');
end


