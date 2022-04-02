function varargout = HM_WideField(varargin)
% HM_WIDEFIELD MATLAB code for HM_WideField.fig
%      HM_WIDEFIELD, by itself, creates a new HM_WIDEFIELD or raises the existing
%      singleton*.
%
%      H = HM_WIDEFIELD returns the handle to a new HM_WIDEFIELD or the handle to
%      the existing singleton*.
%
%      HM_WIDEFIELD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HM_WIDEFIELD.M with the given input arguments.
%
%      HM_WIDEFIELD('Property','Value',...) creates a new HM_WIDEFIELD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HM_WideField_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HM_WideField_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HM_WideField

% Last Modified by GUIDE v2.5 24-Dec-2015 21:31:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HM_WideField_OpeningFcn, ...
    'gui_OutputFcn',  @HM_WideField_OutputFcn, ...
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


% --- Executes just before HM_WideField is made visible.
function HM_WideField_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HM_WideField (see VARARGIN)

% Choose default command line output for HM_WideField
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HM_WideField wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HM_WideField_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.roiList,'String','');
set(handles.roiList,'UserData',[]);

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose File','Multiselect','on');
cd(tiff_path)
if ~iscell(tiff_filename);
    tiff_filename = {tiff_filename};
end

set(handles.tiffFile,'string',[tiff_path tiff_filename{1}]);
% Spawn new image window if new or no match
if isfield(handles,'guiAxes') == 0
    % Set up figure for keyboard shortcuts
    handles.tiff_fig = figure('KeyPressFcn', {@keyPress, hObject}); 
    % Set up the mouse wheel for image scrolling
    set(handles.tiff_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, hObject});
    handles.guiAxes = axes;
    % Create scroll bar
    ypos = [0 0 1 0.05];
    handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
    listener_imgSlider = addlistener(handles.imgSlider,'Action',@imgSlider_Listener);
    set(listener_imgSlider,'CallbackTarget',hObject);
    guidata(hObject,handles);
end
if isfield(handles,'guiAxes') == 1
    if any(ismember(findall(0,'type','axes'),handles.guiAxes)) == 0;
        % Set up figure for keyboard shortcuts
        handles.tiff_fig = figure('KeyPressFcn', {@keyPress, hObject}); 
        % Set up the mouse wheel for image scrolling
        set(handles.tiff_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, hObject});
        handles.guiAxes = axes;
        % Create scroll bar
        ypos = [0 0 1 0.05];
        handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
        listener_imgSlider = addlistener(handles.imgSlider,'Action',@imgSlider_Listener);
        set(listener_imgSlider,'CallbackTarget',hObject);
        guidata(hObject,handles);
    end
end

% Initialize concatenated matrix by finding total number of frames
total_numframes = 0;
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    imageinfo = [];
    imageinfo=imfinfo([img_filename],'tiff');
    numframes=length(imageinfo);
    total_numframes = total_numframes + numframes;
end

M=imageinfo(1).Width;
N=imageinfo(1).Height;
im_concat = zeros(N*M,total_numframes,'uint16');

curr_frame = 1;
for i = 1:length(tiff_filename);
    
    img_filename = [tiff_path tiff_filename{i}];
    
    % Load in image file
    imageinfo = [];
    imageinfo=imfinfo([img_filename],'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    im = [];
    im_temp = [];
    disp('Loading file....')
        for loadframe = 1:numframes
            im_temp = imread([img_filename],'tiff',loadframe,'Info',imageinfo);
            im_temp = im_temp(:);
            im_concat(:,curr_frame) = im_temp;
            curr_frame = curr_frame+1;
            disp(['Frame ' num2str(loadframe) '/' num2str(numframes)]);
        end
end

handles.im = im_concat;

axes(handles.guiAxes);
cla(handles.guiAxes);
imagesc(reshape(handles.im(:,1),N,M),'CDataMapping','scaled');
xlim([0 M]);
ylim([0 N]);
colormap(gray)

caxis([1000 8000]);
set(handles.imgMin,'string',1000);
set(handles.imgMax,'string',8000);
numframes=size(handles.im,2);
set(handles.imgSlider,'Min',1);
if numframes ~= 1;
    set(handles.imgSlider,'Enable','on');
    set(handles.imgSlider,'Max',numframes);
    set(handles.imgSlider,'Value',1);
    set(handles.imgSlider,'SliderStep',[1/(numframes-1), 10/(numframes-1)]);
else
    set(handles.imgSlider,'Enable','off');
end

% Store the image size
handles.N = N;
handles.M = M;

% Make sure figure toolbar is available
if isfield(handles,'tiff_fig')
    set(handles.tiff_fig,'toolbar','figure','Color','w','Position',[300,300,500,500]);
end

% Update handles structure
guidata(hObject, handles);


function tiffFile_Callback(hObject, eventdata, handles)
% hObject    handle to tiffFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tiffFile as text
%        str2double(get(hObject,'String')) returns contents of tiffFile as a double

set(handles.roiList,'String','');
set(handles.roiList,'UserData',[]);

[tiff_fullfilename] = get(handles.tiffFile,'string');
im_s(:,:)=double(imread([tiff_fullfilename],'tiff',1));
axes(handles.guiAxes);
cla(handles.guiAxes);
imagesc(im_s,'CDataMapping','scaled');
colormap(gray)
[imgMin imgMax] = caxis;
set(handles.imgMin,'string',imgMin);
set(handles.imgMax,'string',imgMax);
imageinfo=imfinfo(tiff_fullfilename,'tiff');
numframes=length(imageinfo);
set(handles.imgSlider,'Min',1);
if numframes ~= 1;
    set(handles.imgSlider,'Enable','on');
    set(handles.imgSlider,'Max',numframes);
    set(handles.imgSlider,'Value',1);
    set(handles.imgSlider,'SliderStep',[1/(numframes-1), 10/(numframes-1)]);
else
    set(handles.imgSlider,'Enable','off');
end


% --- Executes during object creation, after setting all properties.
function tiffFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tiffFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in createROI.
function createROI_Callback(hObject, eventdata, handles)
% hObject    handle to createROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
n_polygon = size(get(handles.roiList,'string'),1) + 1;

% Draw ROI polygon
[cellMask,polygon.ROI{n_polygon}(:,1),polygon.ROI{n_polygon}(:,2)] = roipoly;

% Make the polygon closed if it's not already
if polygon.ROI{n_polygon}(1,1) ~= polygon.ROI{n_polygon}(end,1)
    polygon.ROI{n_polygon} = [polygon.ROI{n_polygon};polygon.ROI{n_polygon}(1,:)];
end

% Display just an outline when not selected
polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1), polygon.ROI{n_polygon}(:,2));
set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});

set(handles.roiList,'string',num2str((1:n_polygon)'));
set(handles.roiList,'UserData',polygon);


% --- Executes on button press in deleteROI.
function deleteROI_Callback(hObject, eventdata, handles)
% hObject    handle to deleteROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.guiAxes);
listselect = get(handles.roiList,'value');
polygon = get(handles.roiList,'UserData');
n_polygon = size(get(handles.roiList,'string'),1);

delete(polygon.handles{listselect});
polygon.ROI(listselect) = [];
polygon.handles(listselect) = [];

set(handles.roiList,'value',1);
set(handles.roiList,'string',num2str((1:n_polygon-1)'));
set(handles.roiList,'UserData',polygon);


% --- Executes on button press in saveROI.
function saveROI_Callback(hObject, eventdata, handles)
% hObject    handle to saveROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveROI

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);

[roi_savefile roi_savepath] = uiputfile('.roi','Create save filename',tiff_filename);

imageinfo=imfinfo(tiff_fullfilename,'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

polygon = get(handles.roiList,'UserData');

% Make sure positions are updated in case one is currently selected
axes(handles.guiAxes)
oldpoly = cellfun('isclass', polygon.handles,'impoly');
if any(oldpoly)
    oldROI = find(oldpoly == 1);
    polygon.ROI{oldROI} = getPosition(polygon.handles{oldROI});
    delete(polygon.handles{oldROI})
    polygon.handles{oldROI} = line(polygon.ROI{oldROI}(:,1), polygon.ROI{oldROI}(:,2));
    set(polygon.handles{oldROI}, 'ButtonDownFcn', {@selectROI, handles});
end

% Update polygon handles
set(handles.roiList,'UserData',polygon);

% Don't bother saving handles
polygon = rmfield(polygon,'handles');
if isfield(polygon,'bghandles');
    polygon = rmfield(polygon,'bghandles')
end

save([roi_savepath roi_savefile],'polygon')%'roi_trace'


% --- Executes on button press in loadROI.
function loadROI_Callback(hObject, eventdata, handles)
% hObject    handle to loadROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);

[roi_file roi_path] = uigetfile('.roi','Pick ROI File',tiff_filename);

imageinfo=imfinfo(tiff_fullfilename,'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

axes(handles.guiAxes);

% Delete existing polygons and lines
currpoly = findobj(handles.guiAxes,'tag','impoly');
delete(currpoly)
currpoly = findobj(handles.guiAxes,'type','line');
delete(currpoly)

load([roi_path roi_file],'-MAT')

for n_polygon = 1:length(polygon.ROI)   
    % Make the polygon closed if it's not already
    if polygon.ROI{n_polygon}(1,1) ~= polygon.ROI{n_polygon}(end,1)
        polygon.ROI{n_polygon} = [polygon.ROI{n_polygon};polygon.ROI{n_polygon}(1,:)];
    end
    
    polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1),polygon.ROI{n_polygon}(:,2));
    set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});
end

set(handles.roiList,'Value',1);
set(handles.roiList,'string',num2str((1:n_polygon)'));
set(handles.roiList,'UserData',polygon);


function selectROI(hObject, eventdata, handles)
% Execute when clicking on ROI polygon

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');

oldpoly = cellfun('isclass', polygon.handles,'impoly');
if any(oldpoly)
    oldROI = find(oldpoly == 1);
    polygon.ROI{oldROI} = getPosition(polygon.handles{oldROI});
    % make the polygon closed if it's not already
    if polygon.ROI{oldROI}(1,1) ~= polygon.ROI{oldROI}(end,1)
        polygon.ROI{oldROI} = [polygon.ROI{oldROI};polygon.ROI{oldROI}(1,:)];
    end
    delete(polygon.handles{oldROI})
    polygon.handles{oldROI} = line(polygon.ROI{oldROI}(:,1), polygon.ROI{oldROI}(:,2));
    set(polygon.handles{oldROI}, 'ButtonDownFcn', {@selectROI, handles});
end

% Redraw all polygons in blue
currpoly = findobj(handles.guiAxes,'type','line');

for n_polygon = 1:length(polygon.ROI)
    set(currpoly(n_polygon), 'color', 'b');
end

% Draw selected polygon in red
roiClicked = find([polygon.handles{:}] == hObject);
delete(polygon.handles{roiClicked})
if length(polygon.ROI{roiClicked}) < 20
    polygon.handles{roiClicked} = ...
        impoly(gca,[polygon.ROI{roiClicked}(:,1) polygon.ROI{roiClicked}(:,2)]);
    setColor(polygon.handles{roiClicked},'r')
else
    polygon.handles{roiClicked} = line(polygon.ROI{roiClicked}(:,1), polygon.ROI{roiClicked}(:,2),'color','r');
    set(polygon.handles{roiClicked}, 'ButtonDownFcn', {@selectROI, handles});
end

% Update roi list
set(handles.roiList,'value',roiClicked);
set(handles.roiList,'UserData',polygon);


% --- Executes on selection change in roiList.
function roiList_Callback(hObject, eventdata, handles)
% hObject    handle to roiList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roiList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiList

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
listselect = get(handles.roiList,'value');

selectROI(polygon.handles{listselect},eventdata,handles)


% --- Executes during object creation, after setting all properties.
function roiList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function imgFrame_Callback(hObject, eventdata, handles)
% hObject    handle to imgFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgFrame as text
%        str2double(get(hObject,'String')) returns contents of imgFrame as a double

frame = str2num(get(handles.imgFrame,'String'));
set(handles.imgSlider,'Value',frame);


% --- Executes during object creation, after setting all properties.
function imgFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function imgMin_Callback(hObject, eventdata, handles)
% hObject    handle to imgMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgMin as text
%        str2double(get(hObject,'String')) returns contents of imgMin as a double

imgMin = get(handles.imgMin,'string');
imgMax = get(handles.imgMax,'string');
imgMin = str2num(imgMin);
imgMax = str2num(imgMax);
axes(handles.guiAxes);
caxis([imgMin imgMax]);

% --- Executes during object creation, after setting all properties.
function imgMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function imgMax_Callback(hObject, eventdata, handles)
% hObject    handle to imgMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgMax as text
%        str2double(get(hObject,'String')) returns contents of imgMax as a double

imgMin = get(handles.imgMin,'string');
imgMax = get(handles.imgMax,'string');
imgMin = str2num(imgMin);
imgMax = str2num(imgMax);
axes(handles.guiAxes);
caxis([imgMin imgMax]);

% --- Executes during object creation, after setting all properties.
function imgMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in avgImg.
function avgImg_Callback(hObject, eventdata, handles)
% hObject    handle to avgImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of avgImg

avgImg = get(handles.avgImg,'value');

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

if avgImg == 0;
    set(handles.imgSlider,'Enable','On');
    tiff_fullfilename = get(handles.tiffFile,'string');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    imagesc(reshape(handles.im(:,1),N,M),'CDataMapping','scaled');
    colormap(gray)
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
    
elseif avgImg == 1;
    set(handles.imgSlider,'Enable','Off');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    imagesc(reshape(nanmean(handles.im,2),N,M),'CDataMapping','scaled');
    colormap(gray)
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
end


% --- Executes on button press in getROIPixel.
function getROIPixel_Callback(hObject, eventdata, handles)
% hObject    handle to getROIPixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);

[pixel_savefile pixel_savepath] = uiputfile('.pixel','Create save filename',tiff_filename);

polygon = get(handles.roiList,'UserData');

% Define region of interest
xv = polygon.ROI{1}(:,1)/4; % Downsampled by 4
yv = polygon.ROI{1}(:,2)/4; % Downsampled by 4

% Get 128x128 grid for imaging window
xq = repmat([1:128]',1,128)';
yq = repmat([1:128]',1,128);
[in,on] = inpolygon(xq,yq,xv,yv);

roiPixelInNum = find(reshape(in,128*128,1)); % pixel number inside roi
roiPixelOnNum = find(reshape(on,128*128,1)); % pixel number outside roi
roiPixelNum = union(roiPixelInNum,roiPixelOnNum);

save([pixel_savepath pixel_savefile],'roiPixelNum')


% --- Executes on button press in getCoodinate.
function getCoodinate_Callback(hObject, eventdata, handles)
% hObject    handle to getCoodinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);

[coordinatePixel_savefile coordinatePixel_savepath] = uiputfile('.coordinatePixel','Create save filename',tiff_filename);

set(0,'currentfigure',handles.tiff_fig)
uiwait(msgbox('Click bregma and lambda'));
[x,y] = ginput(2);
bregma = [x(1),y(1)];
lambda = [x(2),y(2)];
coordinate{1} = bregma;
coordinate{2} = lambda;
bregmaToLambdaPixel = ((x(1)-x(2))^2+(y(1)-y(2))^2)^0.5;
bregmaToLambdaLength = 4.21; % (mm)
scaleFac = bregmaToLambdaPixel/bregmaToLambdaLength;
hold on
plot(bregma(1),bregma(2),'o','Color','w','MarkerSize',4);
plot(lambda(1),lambda(2),'o','Color','w','MarkerSize',4);

save([coordinatePixel_savepath coordinatePixel_savefile],'coordinate')

%%% ADD SLIDE FUNCTION
%%% ATOHA GOOD