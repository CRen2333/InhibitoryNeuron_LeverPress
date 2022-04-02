function varargout = GUI_CR_WF_PupilCheck_Fitting(varargin)
% GUI_CR_WF_PUPILCHECK_FITTING MATLAB code for GUI_CR_WF_PupilCheck_Fitting.fig
%      GUI_CR_WF_PUPILCHECK_FITTING, by itself, creates a new GUI_CR_WF_PUPILCHECK_FITTING or raises the isfielding
%      singleton*.
%
%      H = GUI_CR_WF_PUPILCHECK_FITTING returns the handle to a new GUI_CR_WF_PUPILCHECK_FITTING or the handle to
%      the isfielding singleton*.
%
%      GUI_CR_WF_PUPILCHECK_FITTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CR_WF_PUPILCHECK_FITTING.M with the given input arguments.
%
%      GUI_CR_WF_PUPILCHECK_FITTING('Property','Value',...) creates a new GUI_CR_WF_PUPILCHECK_FITTING or raises the
%      isfielding singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_CR_WF_PupilCheck_Fitting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_CR_WF_PupilCheck_Fitting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_CR_WF_PupilCheck_Fitting

% Last Modified by GUIDE v2.5 06-Nov-2017 11:18:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_CR_WF_PupilCheck_Fitting_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_CR_WF_PupilCheck_Fitting_OutputFcn, ...
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


% --- Executes just before GUI_CR_WF_PupilCheck_Fitting is made visible.
function GUI_CR_WF_PupilCheck_Fitting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_CR_WF_PupilCheck_Fitting (see VARARGIN)

% Choose default command line output for GUI_CR_WF_PupilCheck_Fitting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_CR_WF_PupilCheck_Fitting wait for user response (see UIRESUME)
% uiwait(handles.fig);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_CR_WF_PupilCheck_Fitting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_Load.
function pushbutton_Load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Initialize

[handles.FileName,handles.PathName] = uigetfile('*.mat','Select the Movie file');
load([handles.PathName handles.FileName]);
handles.Image = Eye_ROI_filtered;
handles.bestFits = bestFits;
handles.bestFits_refine = bestFits_refine;
handles.PupilDia_refine = PupilDia_refine;
handles.Filename = [Initial '_' Date '_' Animal '_Fitting'];
set(handles.edit_path,'string',[handles.FileName]);
set(handles.edit_TotalFrameNum,'string',num2str(size(Eye_ROI_filtered,3)));
set(handles.edit_FrameRate,'String','30');
set(handles.edit_Frame,'string','1');

handles.TotalFrameSelected = nan(1,size(Eye_ROI_filtered,3)); % default: nan = defult fitting;
set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));

handles.FrameNum = str2double(get(handles.edit_Frame,'String'));
axes(handles.axes_Image);
cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal; axis off;
hold on;
% handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
% set(handles.h1,'color',[0.75 0.93 1]);
% handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
% set(handles.h2,'color',[0.7 0.9 0.25]);
% handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
% set(handles.h3,'color',[0.93 0.69 0.13]);
handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
set(handles.h4,'color',[1 0.7 0.7]);
            
set(handles.slider_Image,'Min',1);
set(handles.slider_Image,'Enable','on');
set(handles.slider_Image,'Value',1);
set(handles.slider_Image,'Max',size(handles.Image,3));
set(handles.slider_Image,'SliderStep',[1/(size(handles.Image,3)-1) 10/(size(handles.Image,3)-1)]);

set(handles.checkbox_Fit_1,'Value',0);
set(handles.checkbox_Fit_2,'Value',0);
set(handles.checkbox_Fit_3,'Value',0);
set(handles.checkbox_Bad_Fitting,'Value',0);

axes(handles.axes_PupilDia);
cla; set(handles.axes_PupilDia,'color',[0.31 0.31 0.31]);
plot(handles.PupilDia_refine);
hold on; line([handles.FrameNum handles.FrameNum],ylim,'color','w');

assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
guidata(hObject,handles);


function edit_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path as text
%        str2double(get(hObject,'String')) returns contents of edit_path as a double
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_Image_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.FrameNum = round(get(handles.slider_Image,'Value'));
set(handles.edit_Frame,'string',num2str(handles.FrameNum));
axes(handles.axes_PupilDia);
cla; set(handles.axes_PupilDia,'color',[0.31 0.31 0.31]);
plot(handles.PupilDia_refine);
hold on; line([handles.FrameNum handles.FrameNum],ylim,'color','w');
if isnan(handles.TotalFrameSelected(handles.FrameNum))
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal; axis off;
    hold on;
%     handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
%     set(handles.h1,'color',[0.75 0.93 1]);
%     handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
%     set(handles.h2,'color',[0.7 0.9 0.25]);
%     handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
%     set(handles.h3,'color',[0.93 0.69 0.13]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 1
    set(handles.checkbox_Fit_1,'Value',1);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
    set(handles.h1,'color',[0.75 0.93 1]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 2
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',1);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
    set(handles.h2,'color',[0.7 0.9 0.25]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 3
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',1);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
    set(handles.h3,'color',[0.93 0.69 0.13]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 0
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',1);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
end
set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));    
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_Image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function edit_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Frame as text
%        str2double(get(hObject,'String')) returns contents of edit_Frame as a double
set(handles.edit_Frame,'string',num2str(handles.FrameNum));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_FrameRate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Frame as text
%        str2double(get(hObject,'String')) returns contents of edit_Frame as a double
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_FrameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_TotalFrameNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_TotalFrameNum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Frame as text
%        str2double(get(hObject,'String')) returns contents of edit_Frame as a double
guidata(hObject,handles);

% --- Executes on button press in checkbox_Fit_1.
function checkbox_Fit_1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Fit_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Fit_1
handles.Fit_1 = get(handles.checkbox_Fit_1,'Value');
if handles.Fit_1 == 1
    cla
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = 1;
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
    set(handles.h1,'color',[0.75 0.93 1]);
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.Fit_1 == 0
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    set(handles.checkbox_Fit_1,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = nan;
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
    rmfield(handles,'h1');
end
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
guidata(hObject,handles);


% --- Executes on button press in checkbox_Fit_2.
function checkbox_Fit_2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Fit_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Fit_2
handles.Fit_2 = get(handles.checkbox_Fit_2,'Value');
if handles.Fit_2 == 1
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = 2;
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
    set(handles.h2,'color',[0.7 0.9 0.25]);
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.Fit_2 == 0
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    set(handles.checkbox_Fit_2,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = nan;
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
end
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
guidata(hObject,handles);

% --- Executes on button press in checkbox_Fit_3.
function checkbox_Fit_3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Fit_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Fit_3
handles.Fit_3 = get(handles.checkbox_Fit_3,'Value');
if handles.Fit_3 == 1
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = 3;
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
    set(handles.h3,'color',[0.93 0.69 0.13]);
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.Fit_3 == 0
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    set(handles.checkbox_Fit_3,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = nan;
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
end
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
guidata(hObject,handles);


function edit_TotalFrameSelected_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TotalFrameSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TotalFrameSelected as text
%        str2double(get(hObject,'String')) returns contents of edit_TotalFrameSelected as a double
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_TotalFrameSelected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TotalFrameSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_Finish.
function pushbutton_Finish_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
TotalFrameSelected = handles.TotalFrameSelected;
save([handles.PathName handles.FileName],'TotalFrameSelected','-append');
h = msgbox('Saving Completed');
guidata(hObject,handles);


% --- Executes on key press with focus on fig or any of its controls.
function fig_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fig (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
frame = handles.FrameNum;
numframes = size(handles.Image,3);
switch eventdata.Key
    case 'leftarrow'
        frame = round(frame)-1;       
    case 'rightarrow'
        frame = round(frame)+1;
end
% constrain scroll movement to number of frames
if frame < 1
    frame = 1;
elseif frame > numframes
    frame = numframes;
end
handles.FrameNum = frame;
set(handles.slider_Image,'Value',frame);
set(handles.edit_Frame,'string',num2str(handles.FrameNum));
axes(handles.axes_PupilDia);
cla; set(handles.axes_PupilDia,'color',[0.31 0.31 0.31]);
plot(handles.PupilDia_refine);
hold on; line([handles.FrameNum handles.FrameNum],ylim,'color','w');
if isnan(handles.TotalFrameSelected(handles.FrameNum))
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
%     handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
%     set(handles.h1,'color',[0.75 0.93 1]);
%     handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
%     set(handles.h2,'color',[0.7 0.9 0.25]);
%     handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
%     set(handles.h3,'color',[0.93 0.69 0.13]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 1
    set(handles.checkbox_Fit_1,'Value',1);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
    set(handles.h1,'color',[0.75 0.93 1]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 2
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',1);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
    set(handles.h2,'color',[0.7 0.9 0.25]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 3
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',1);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
    set(handles.h3,'color',[0.93 0.69 0.13]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 0
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',1);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
end
set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));    
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);

guidata(hObject,handles);


% --- Executes on scroll wheel click while the figure is in focus.
function fig_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to fig (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
frame = handles.FrameNum;
frame_move = eventdata.VerticalScrollCount;
numframes = size(handles.Image,3);
frame = frame+frame_move;
% constrain scroll movement to number of frames
if frame < 1
    frame = 1;
elseif frame > numframes
    frame = numframes;
end
handles.FrameNum = frame;
set(handles.slider_Image,'Value',frame);
set(handles.edit_Frame,'string',num2str(handles.FrameNum));
axes(handles.axes_PupilDia);
cla; set(handles.axes_PupilDia,'color',[0.31 0.31 0.31]);
plot(handles.PupilDia_refine);
hold on; line([handles.FrameNum handles.FrameNum],ylim,'color','w');
if isnan(handles.TotalFrameSelected(handles.FrameNum))
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
%     handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
%     set(handles.h1,'color',[0.75 0.93 1]);
%     handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
%     set(handles.h2,'color',[0.7 0.9 0.25]);
%     handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
%     set(handles.h3,'color',[0.93 0.69 0.13]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 1
    set(handles.checkbox_Fit_1,'Value',1);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h1 = ellipse(handles.bestFits{handles.FrameNum}(1,3),handles.bestFits{handles.FrameNum}(1,4),handles.bestFits{handles.FrameNum}(1,5)*pi/180,handles.bestFits{handles.FrameNum}(1,1),handles.bestFits{handles.FrameNum}(1,2));
    set(handles.h1,'color',[0.75 0.93 1]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 2
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',1);
    set(handles.checkbox_Fit_3,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h2 = ellipse(handles.bestFits{handles.FrameNum}(2,3),handles.bestFits{handles.FrameNum}(2,4),handles.bestFits{handles.FrameNum}(2,5)*pi/180,handles.bestFits{handles.FrameNum}(2,1),handles.bestFits{handles.FrameNum}(2,2));
    set(handles.h2,'color',[0.7 0.9 0.25]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 3
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',1);
    set(handles.checkbox_Bad_Fitting,'Value',0);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    hold on;
    handles.h3 = ellipse(handles.bestFits{handles.FrameNum}(3,3),handles.bestFits{handles.FrameNum}(3,4),handles.bestFits{handles.FrameNum}(3,5)*pi/180,handles.bestFits{handles.FrameNum}(3,1),handles.bestFits{handles.FrameNum}(3,2));
    set(handles.h3,'color',[0.93 0.69 0.13]);
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
elseif handles.TotalFrameSelected(handles.FrameNum) == 0
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    set(handles.checkbox_Bad_Fitting,'Value',1);
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
end
set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));    
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);

guidata(hObject,handles);


% --- Executes on button press in checkbox_Bad_Fitting.
function checkbox_Bad_Fitting_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Bad_Fitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Bad_Fitting
handles.Bad_Fitting = get(handles.checkbox_Bad_Fitting,'Value');
if handles.Bad_Fitting == 1
    set(handles.checkbox_Fit_1,'Value',0);
    set(handles.checkbox_Fit_2,'Value',0);
    set(handles.checkbox_Fit_3,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = 0;
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
elseif handles.Bad_Fitting == 0
    set(handles.checkbox_Bad_Fitting,'Value',0);
    handles.TotalFrameSelected(handles.FrameNum) = nan;
    set(handles.edit_TotalFrameSelected,'String',num2str(sum(~isnan(handles.TotalFrameSelected))));
    axes(handles.axes_Image);
    cla;imagesc(handles.Image(:,:,handles.FrameNum)/255, [0,1]); colormap(gray); axis equal;
    handles.h4 = ellipse(handles.bestFits_refine{handles.FrameNum}(1,3),handles.bestFits_refine{handles.FrameNum}(1,4),handles.bestFits_refine{handles.FrameNum}(1,5)*pi/180,handles.bestFits_refine{handles.FrameNum}(1,1),handles.bestFits_refine{handles.FrameNum}(1,2));
    set(handles.h4,'color',[1 0.7 0.7]);
end
assignin('base','TotalFrameSelected',handles.TotalFrameSelected);
guidata(hObject,handles);
