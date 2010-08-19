function varargout = prColormapPanel(varargin)
% PRCOLORMAPPANEL M-file for prColormapPanel.fig
%      PRCOLORMAPPANEL, by itself, creates a new PRCOLORMAPPANEL or raises the existing
%      singleton*.
%
%      H = PRCOLORMAPPANEL returns the handle to a new PRCOLORMAPPANEL or the handle to
%      the existing singleton*.
%
%      PRCOLORMAPPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRCOLORMAPPANEL.M with the given input arguments.
%
%      PRCOLORMAPPANEL('Property','Value',...) creates a new PRCOLORMAPPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prColormapPanel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prColormapPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help prColormapPanel

% Last Modified by GUIDE v2.5 23-Feb-2005 20:40:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prColormapPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @prColormapPanel_OutputFcn, ...
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


% --- Executes just before prColormapPanel is made visible.
function prColormapPanel_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% check if the colormap was set previsouly
hfsmC=findall(0,'Tag','figure1','Name','Cell dynamics analysis');
if ~isempty(hfsmC)
    % Get current project from prAlphaPanel
    handlesprAlphaPanel = guidata(hfsmC);
    if isfield(handlesprAlphaPanel,'color_map_parameters')
        color_map_parameters = handlesprAlphaPanel.color_map_parameters;
        set(handles.protrusion,'String',handlesprAlphaPanel.color_map_parameters.protrusion);
        set(handles.poly,'String',handlesprAlphaPanel.color_map_parameters.poly);
        set(handles.protrusion_u,'String',handlesprAlphaPanel.color_map_parameters.protrusion_u);
        set(handles.poly_u,'String',handlesprAlphaPanel.color_map_parameters.poly_u);        
        set(handles.activity1,'String',handlesprAlphaPanel.color_map_parameters.activity1);
        set(handles.activity2,'String',handlesprAlphaPanel.color_map_parameters.activity2);
        set(handles.activity1_u,'String',handlesprAlphaPanel.color_map_parameters.activity1_u);
        set(handles.activity2_u,'String',handlesprAlphaPanel.color_map_parameters.activity2_u);       
        set(handles.de_protrusion,'String',handlesprAlphaPanel.color_map_parameters.de_protrusion);
        set(handles.de_poly,'String',handlesprAlphaPanel.color_map_parameters.de_poly);
        set(handles.de_activity1,'String',handlesprAlphaPanel.color_map_parameters.de_activity1);
        set(handles.de_activity2,'String',handlesprAlphaPanel.color_map_parameters.de_activity2);       
    else
        set(handles.protrusion,'String',-20);
        set(handles.poly,'String',-4e-4);
        set(handles.protrusion_u,'String',20);
        set(handles.poly_u,'String',4e-4);       
        set(handles.activity1,'String',1000);
        set(handles.activity2,'String',1000);
        set(handles.activity1_u,'String',2500);
        set(handles.activity2_u,'String',2500);        
        set(handles.de_protrusion,'String',4);
        set(handles.de_poly,'String',3e-4);
        set(handles.de_activity1,'String',0);
        set(handles.de_activity2,'String',0);
    end
end

%fill the predefined list list
vars{1} = 'Generic';
vars{2} = 'RacPAK 1 (control)';
vars{3} = 'RacPAK 2'; 
vars{4} = 'RacPAK 3'; 
vars{5} = 'Newt 1 (2135)'; 
vars{6} = 'Newt 2 '; 
vars{7} = 'Newt 3 ';
vars{8} = 'TM: s570 ';
vars{9} = 'Cdc42 1';
vars{10} = 'Cdc42 2'; 
vars{11} = 'Cdc42 wound'; 
vars{12} = 'Rho 1 new cell'; 
vars{13} = 'Rho 2'; 
vars{14} = 'Paxilin cell3'; 
vars{15} = 'Paxilin cell4'; 
vars{16} = 'Intermediate filaments 0421';

set(handles.colormap_list,'String',vars);

set(handles.protrusion_u,'Enable','off');
set(handles.poly_u,'Enable','off');
    
guidata(hObject, handles);

set(handles.prcolormapFigure,'Visible','on');
waitfor(handles.prcolormapFigure,'Visible');
% UIWAIT makes prColormapPanel wait for user response (see UIRESUME)
% uiwait(handles.prcolormapFigure);



function varargout = prColormapPanel_OutputFcn(hObject, eventdata, handles) 
color_map_parameters.protrusion = str2double(get(handles.protrusion, 'String'));
color_map_parameters.poly       = str2double(get(handles.poly, 'String'));
color_map_parameters.activity1  = str2double(get(handles.activity1, 'String'));
color_map_parameters.activity2  = str2double(get(handles.activity2, 'String'));

color_map_parameters.protrusion_u = str2double(get(handles.protrusion_u, 'String'));
color_map_parameters.poly_u       = str2double(get(handles.poly_u, 'String'));
color_map_parameters.activity1_u  = str2double(get(handles.activity1_u, 'String'));
color_map_parameters.activity2_u  = str2double(get(handles.activity2_u, 'String'));

color_map_parameters.de_protrusion  = str2double(get(handles.de_protrusion, 'String'));
color_map_parameters.de_poly        = str2double(get(handles.de_poly, 'String'));
color_map_parameters.de_activity1   = str2double(get(handles.de_activity1, 'String'));
color_map_parameters.de_activity2   = str2double(get(handles.de_activity2, 'String'));
varargout{1}    =   color_map_parameters;
delete(hObject);


function colormap_list_Callback(hObject, eventdata, handles)
% vars{1} = 'Generic';
% vars{2} = 'RacPAK 1 (control)';
% vars{3} = 'RacPAK 2'; 
% vars{4} = 'RacPAK 3'; 
% vars{5} = 'Newt 1 (2135)'; 
% vars{6} = 'Newt 2 '; 
% vars{7} = 'Newt 3 '; 
% vars{8} = 'TM: 570'; 
% vars{9} = 'Cdc42 1';
% vars{10} = 'Cdc42 2'; 
% vars{11} = 'Cdc42 wound'; 
% vars{12} = 'Rho 1 new cell'; 
% vars{13} = 'Rho 2'; 
% vars{14} = 'Paxilin cell3'; 
% vars{15} = 'Paxilin cell4'; 
% vars{16} = 'Intermediate filaments 0421';


val = get(hObject,'Value');
if val < 9
    set(handles.activity1,'String',0);
    set(handles.activity2,'String',0);
    set(handles.activity1_u,'String',0);
    set(handles.activity2_u,'String',0);    
    set(handles.de_activity1,'String',0);
    set(handles.de_activity2,'String',0);     
else
    set(handles.poly,'String',0);  
    set(handles.de_poly,'String',0);
end

if val == 1
    set(handles.protrusion,'String',-15);
    set(handles.poly,'String',-3e-4);
    set(handles.protrusion_u,'String',15);
    set(handles.poly_u,'String',3e-4);    
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4);  
elseif val == 2
    set(handles.protrusion,'String',-25);
    set(handles.poly,'String',-4e-4);
    set(handles.protrusion_u,'String',25);
    set(handles.poly_u,'String',4e-4);    
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4);    
elseif val == 3
    set(handles.protrusion,'String',-20);
    set(handles.poly,'String',-4e-4);
    set(handles.protrusion_u,'String',20);
    set(handles.poly_u,'String',4e-4);    
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4); 
elseif val == 4
    set(handles.protrusion,'String',-20);
    set(handles.poly,'String',-4e-4);  
    set(handles.protrusion_u,'String',20);
    set(handles.poly_u,'String',4e-4);      
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4);  
elseif val == 5
    set(handles.protrusion,'String',-20);
    set(handles.protrusion_u,'String',20);    
    set(handles.poly,'String',-4e-4);
    set(handles.poly_u,'String',4e-4);    
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4);       
elseif val == 6
    set(handles.protrusion,'String',-20);
    set(handles.protrusion_u,'String',20);    
    set(handles.poly,'String',-4e-4);
    set(handles.poly_u,'String',4e-4);    
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4);   
 elseif val == 7
    set(handles.protrusion,'String',-20);
    set(handles.protrusion_u,'String',20);    
    set(handles.poly,'String',-4e-4);
    set(handles.poly_u,'String',4e-4);    
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4);      
 elseif val == 8
    set(handles.protrusion,'String',-80);
    set(handles.protrusion_u,'String',80);    
    set(handles.poly,'String',-6e-4);
    set(handles.poly_u,'String',6e-4);    
    set(handles.de_protrusion,'String',4);
    set(handles.de_poly,'String',3e-4);   
elseif val == 9
    set(handles.protrusion,'String',-15);
    set(handles.activity1,'String',1000);
    set(handles.activity2,'String',1000);   
    set(handles.activity1_u,'String',3500);
    set(handles.activity2_u,'String',3500);        
    set(handles.de_protrusion,'String',4);
    set(handles.de_activity1,'String',2);
    set(handles.de_activity2,'String',2);     
elseif val == 10
    set(handles.protrusion,'String',-20);
    set(handles.activity1,'String',1000);
    set(handles.activity2,'String',1000);   
    set(handles.activity1_u,'String',3500);
    set(handles.activity2_u,'String',3500);      
    set(handles.de_protrusion,'String',4);
    set(handles.de_activity1,'String',3);
    set(handles.de_activity2,'String',3);     
elseif val == 11
    set(handles.protrusion,'String',-20);
    set(handles.activity1,'String',1500);
    set(handles.activity2,'String',1500);  
    set(handles.activity1_u,'String',3500);
    set(handles.activity2_u,'String',3500);      
    set(handles.de_protrusion,'String',20);
    set(handles.de_activity1,'String',100);
    set(handles.de_activity2,'String',4);     
elseif val == 12
    set(handles.protrusion,'String',-20);
    set(handles.activity1,'String',1200);
    set(handles.activity2,'String',1000);   
    set(handles.activity1_u,'String',2600);
    set(handles.activity2_u,'String',3500);     
    set(handles.de_protrusion,'String',4);
    set(handles.de_activity1,'String',2);
    set(handles.de_activity2,'String',2);     
elseif val == 13
    set(handles.protrusion,'String',-20);
    set(handles.activity1,'String',0);
    set(handles.activity2,'String',0);    
    set(handles.activity1_u,'String',3500);
    set(handles.activity2_u,'String',3500);      
    set(handles.de_protrusion,'String',4);
    set(handles.de_activity1,'String',2);
    set(handles.de_activity2,'String',2);   
 elseif val == 14
    set(handles.protrusion,'String',-15);
    set(handles.protrusion_u,'String',15);
    set(handles.poly,'String',-4e-4);
    set(handles.poly_u,'String',4e-4);    
    set(handles.activity1,'String',445);
    set(handles.activity2,'String',445);    
    set(handles.activity1_u,'String',455);
    set(handles.activity2_u,'String',455);      
    set(handles.de_protrusion,'String',4);
    set(handles.de_activity1,'String',2);
    set(handles.de_activity2,'String',2);   
 elseif val == 15
    set(handles.protrusion,'String',-7);
    set(handles.protrusion_u,'String',7);
    set(handles.poly,'String',-4e-4);
    set(handles.poly_u,'String',4e-4);        
    set(handles.activity1,'String',445);
    set(handles.activity2,'String',445);    
    set(handles.activity1_u,'String',457);
    set(handles.activity2_u,'String',457);      
    set(handles.de_protrusion,'String',4);
    set(handles.de_activity1,'String',2);
    set(handles.de_activity2,'String',2);  
 elseif val == 16
    set(handles.protrusion,'String',-15);
    set(handles.protrusion_u,'String',15);
    set(handles.poly,'String',-4e-4);
    set(handles.poly_u,'String',4e-4);        
    set(handles.activity1,'String',60);
    set(handles.activity2,'String',60);    
    set(handles.activity1_u,'String',115);
    set(handles.activity2_u,'String',115);      
    set(handles.de_protrusion,'String',4);
    set(handles.de_activity1,'String',2);
    set(handles.de_activity2,'String',2);    
end

function colormap_list_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function protrusion_Callback(hObject, eventdata, handles)
value = -str2double(get(handles.protrusion,'String'));
set(handles.protrusion_u,'String',num2str(value));
guidata(hObject, handles);

function protrusion_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function poly_Callback(hObject, eventdata, handles)
value = -str2double(get(handles.poly,'String'));
set(handles.poly_u,'String',num2str(value));
guidata(hObject, handles);


function poly_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function flow_Callback(hObject, eventdata, handles)

function flow_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function activity1_Callback(hObject, eventdata, handles)

function activity1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function activity2_Callback(hObject, eventdata, handles)

function activity2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function done_Callback(hObject, eventdata, handles)
prColormapPanel_CloseRequestFcn(hObject, eventdata, handles);

function prColormapPanel_CloseRequestFcn(hObject, eventdata, handles)
% hide the figure here, it will be deleted in the OutputFcn
set(handles.prcolormapFigure,'Visible','off');



function de_protrusion_Callback(hObject, eventdata, handles)
function de_protrusion_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function de_poly_Callback(hObject, eventdata, handles)
function de_poly_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function de_activity2_Callback(hObject, eventdata, handles)
function de_activity2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function protrusion_u_Callback(hObject, eventdata, handles)



function protrusion_u_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function poly_u_Callback(hObject, eventdata, handles)
function poly_u_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function activity1_u_Callback(hObject, eventdata, handles)
function activity1_u_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function activity2_u_Callback(hObject, eventdata, handles)
function activity2_u_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


