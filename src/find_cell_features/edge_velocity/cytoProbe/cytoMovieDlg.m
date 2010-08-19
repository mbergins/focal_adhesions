function varargout = cytoMovieDlg(varargin)
% CYTOMOVIEDLG M-file for cytoMovieDlg.fig
%      CYTOMOVIEDLG by itself, creates a new CYTOMOVIEDLG or raises the
%      existing singleton*.
%
%      H = CYTOMOVIEDLG returns the handle to a new CYTOMOVIEDLG or the handle to
%      the existing singleton*.
%
%      CYTOMOVIEDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYTOMOVIEDLG.M with the given input arguments.
%
%      CYTOMOVIEDLG('Property','Value',...) creates a new CYTOMOVIEDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cytoMovieDlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cytoMovieDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cytoMovieDlg

% Last Modified by GUIDE v2.5 12-Oct-2006 14:27:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cytoMovieDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @cytoMovieDlg_OutputFcn, ...
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

% --- Executes just before cytoMovieDlg is made visible.
function cytoMovieDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cytoMovieDlg (see VARARGIN)


% set the colorbar options
color_map = cell(1,1);
color_map(1) = cellstr('gray');
color_map(2) = cellstr('hot');
color_map(3) = cellstr('cool');
color_map(4) = cellstr('summer');
color_map(5) = cellstr('jet');
color_map(6) = cellstr('hsv');
set(handles.color_map   ,'String', color_map);

bit_depth = cell(1,1);
bit_depth(1) = cellstr('8');
bit_depth(2) = cellstr('10');
bit_depth(3) = cellstr('12');
bit_depth(4) = cellstr('14');
bit_depth(5) = cellstr('16');
set(handles.bit_depth   ,'String', bit_depth);

% set the box color options
box_color = cell(1,1);
box_color(1) = cellstr('white');
box_color(2) = cellstr('black');
box_color(3) = cellstr('red');
box_color(4) = cellstr('green');
box_color(5) = cellstr('blue');
set(handles.box_color   ,'String', box_color);


% Choose default command line output for cytoMovieDlg
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

% UIWAIT makes cytoMovieDlg wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = cytoMovieDlg_OutputFcn(hObject, eventdata, handles)
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
guidata(hObject, handles);
val  = get(handles.color_map,'Value');
list = get(handles.color_map,'String');
color_map = char(list{val}); 

val  = get(handles.bit_depth,'Value');
list = get(handles.bit_depth,'String');
bit_depth = char(list{val}); 

fps  = get(handles.fps,'String');
offset  = get(handles.offset,'String');
scaling = get(handles.scaling,'String');


box  = get(handles.box,'Value');
val  = get(handles.box_color,'Value');
list = get(handles.box_color,'String');
box_color = char(list{val}); 

bar         = get(handles.scale_bar,'Value');
bar_length  = get(handles.bar_length,'String');

timer           = get(handles.timer,'Value');
timer_interval  = get(handles.interval,'String');

automatic_adjust = get(handles.automatic_adjust,'Value');

img_resize = get(handles.img_resize,'String');

handles.output.fps          = str2double(fps);
handles.output.scaling      = str2double(scaling);
handles.output.offset       = str2double(offset);
handles.output.colormap     = color_map;
handles.output.bit_depth    = str2num(bit_depth);
handles.output.box          = box;
handles.output.box_color    = box_color;
handles.output.bar          = bar;
handles.output.bar_length   = str2num(bar_length);
handles.output.automatic_adjust = automatic_adjust;
handles.output.img_resize   = str2num(img_resize);
handles.output.timer            = timer;
handles.output.timer_interval   = str2num(timer_interval);

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


function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function background_Callback(hObject, eventdata, handles)
function roi_Callback(hObject, eventdata, handles)
function timer_Callback(hObject, eventdata, handles)
function nr_user_bands_Callback(hObject, eventdata, handles)
function nr_user_bands_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pp_Callback(hObject, eventdata, handles)
if get(handles.pp, 'Value');
   set(handles.normalise_end, 'Value',0);
end

function normalise_end_Callback(hObject, eventdata, handles)
if get(handles.normalise_end, 'Value');
    set(handles.pp, 'Value',0);
end

function color_map_Callback(hObject, eventdata, handles)
function color_map_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scaling_Callback(hObject, eventdata, handles)
function scaling_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scale_bar_Callback(hObject, eventdata, handles)

function bit_depth_Callback(hObject, eventdata, handles)
function bit_depth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fps_Callback(hObject, eventdata, handles)
function fps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function box_Callback(hObject, eventdata, handles)
function box_color_Callback(hObject, eventdata, handles)

function box_color_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bar_length_Callback(hObject, eventdata, handles)
function bar_length_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function automatic_adjust_Callback(hObject, eventdata, handles)
if get(handles.automatic_adjust,'Value')
    set(handles.offset,'Enable','off');
    set(handles.scaling,'Enable','off');
else
    set(handles.offset,'Enable','on');
    set(handles.scaling,'Enable','on');    
end


function img_resize_Callback(hObject, eventdata, handles)

function img_resize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function interval_Callback(hObject, eventdata, handles)
% hObject    handle to interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interval as text
%        str2double(get(hObject,'String')) returns contents of interval as a double


% --- Executes during object creation, after setting all properties.
function interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


