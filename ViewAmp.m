function varargout = ViewAmp(varargin)
%VIEWAMP M-file for ViewAmp.fig
%      VIEWAMP, by itself, creates a new VIEWAMP or raises the existing
%      singleton*.
%
%      H = VIEWAMP returns the handle to a new VIEWAMP or the handle to
%      the existing singleton*.
%
%      VIEWAMP('Property','Value',...) creates a new VIEWAMP using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ViewAmp_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      VIEWAMP('CALLBACK') and VIEWAMP('CALLBACK',hObject,...) call the
%      local function named CALLBACK in VIEWAMP.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewAmp

% Last Modified by GUIDE v2.5 15-Nov-2019 10:59:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ViewAmp_OpeningFcn, ...
    'gui_OutputFcn',  @ViewAmp_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before ViewAmp is made visible.
function ViewAmp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ViewAmp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ViewAmp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ViewAmp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4


% --- Executes on slider movement.
function ImageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% Figure out which image to show
index = ceil(get(hObject, 'Value'));

% Update existing image object in the GUI using this image data
set(handles.image, 'CData', handles.intensity{index});


% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function TimeInterval_Callback(hObject, eventdata, handles)
% hObject    handle to TimeInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeInterval as text
%        str2double(get(hObject,'String')) returns contents of TimeInterval as a double

% get 'Time interval (s)'
TimeInterval = str2double(get(hObject, 'String'));
handles.TimeInterval = TimeInterval;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function TimeInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SampleRate as text
%        str2double(get(hObject,'String')) returns contents of SampleRate as a double

% get 'Sample rate (fps)'
Fs = str2double(get(hObject, 'String'));
handles.Fs = Fs;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function SampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EquipNiUSB.
function EquipNiUSB_Callback(hObject, eventdata, handles)
% hObject    handle to EquipNiUSB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


s = daq.createSession('ni');
[ch1, idx1] = addAnalogInputChannel(s, 'Dev1', 0, 'Voltage');
[ch2, idx2] = addAnalogInputChannel(s, 'Dev1', 1, 'Voltage');
[ch3, idx3] = addAnalogInputChannel(s, 'Dev1', 7, 'Voltage');
handles.s = s;
guidata(hObject, handles);

warndlg('Ni-Daq is ready!')

% --- Executes on button press in StartTimer.
function StartTimer_Callback(hObject, eventdata, handles)
% hObject    handle to StartTimer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Start timer
s = handles.s;
s.Rate = 1000;
s.DurationInSeconds = handles.NiDaqTime;
% lh = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps, event.Data));
data = startForeground(s);
TimeStamps = (1:size(data, 1))/s.Rate;
guidata(hObject, handles);
% delete(lh);
figure('color', 'w');
plot(TimeStamps, data);

[folder_structure, current_folder] = fileparts(handles.DirectoryName);
if length(data) > 0
    mkdir([folder_structure '\MAT']);
    mkdir([folder_structure '\MAT\' current_folder  ]);
end

SavePath = [folder_structure '\MAT\' handles.expName '_timer.mat'];
save(SavePath, 'data', 'TimeStamps');



function AcqTime_Callback(hObject, eventdata, handles)
% hObject    handle to AcqTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AcqTime as text
%        str2double(get(hObject,'String')) returns contents of AcqTime as a double

% get 'Ni Acquisition time (s)'
NiDaqTime = str2double(get(hObject, 'String'));
handles.NiDaqTime = NiDaqTime;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function AcqTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AcqTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RawPath.
function RawPath_Callback(hObject, eventdata, handles)
% hObject    handle to RawPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DirectoryName = uigetdir();
handles.DirectoryName = DirectoryName;
guidata(hObject,handles);


function PreviewNum_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PreviewNum as text
%        str2double(get(hObject,'String')) returns contents of PreviewNum as a double

% Preview Number
PreviewNum = str2double(get(hObject, 'String'));
handles.PreviewNum = PreviewNum;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PreviewNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PreviewNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load.
function Load_Callback(hObject, ~, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load raws

BeginPoint = 1;
PreviewList = dir([handles.DirectoryName '\*.raw']);

intensity = cell(handles.PreviewNum + 1, 1);
eight_bit = 0; % select 8 bit  or 16 bit raw files; default is 16 bit

% The first image 'intensity0'
fid = fopen([handles.DirectoryName '\' PreviewList(1).name]);
A = fread(fid, 'uint8=>uint8');
fclose(fid);
% allign bits
E = double(A(1:2:end));
F = double(A(2:2:end));
G = 64*E+F/4;
if eight_bit == 1
    intensity0 = reshape(E, [640 480]);
elseif eight_bit == 0
    intensity0 = reshape(G, [640 480])';
end
handles.intensity0 = intensity0;

FileList = dir([handles.DirectoryName '\*.raw']);
if length(FileList) <= handles.PreviewNum
    pause(handles.TimeInterval);
end

for file = BeginPoint:(handles.PreviewNum + BeginPoint)
    
    % read files
    fid = fopen([handles.DirectoryName '\' PreviewList(file).name]);
    A = fread(fid, 'uint8=>uint8');
    fclose(fid);
    % allign bits
    E = double(A(1:2:end));
    F = double(A(2:2:end));
    G = 64*E+F/4;
    if eight_bit == 1
        intensity{file-BeginPoint+1, 1} = reshape(E, [640 480]);
    elseif eight_bit == 0
        intensity{file-BeginPoint+1, 1} = reshape(G, [640 480])';
    end
    
    % Write Files
    %     imwrite(uint16(intensity), [folder_structure '\TIFF\' current_folder '\' PreviewList(file).name '.tiff'],...
    %         'Compression', 'none');
    intensity{file-BeginPoint+1, 1} = intensity{file-BeginPoint+1, 1} - intensity0;
    
end
handles.intensity = intensity;

% Display the first one and store the graphics handle to the imshow object
handles.image = imshow(handles.intensity{1}, 'Parent', handles.axes1);
handles.hZoom = [];

set(handles.ImageSlider, 'Min', 1, 'Max', (handles.PreviewNum + 1), ...
    'SliderStep', [1 1]/(handles.PreviewNum+1 - 1), 'Value', 1)

handles.BeginPoint = BeginPoint+handles.PreviewNum;
guidata(hObject, handles);


% --- Executes on button press in LoadMore.
function LoadMore_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

BeginPoint = handles.BeginPoint;
intensity0 = handles.intensity0;
PreviewList = dir([handles.DirectoryName '\*.raw']);

intensity = cell(handles.PreviewNum + 1, 1);
eight_bit = 0; % select 8 bit  or 16 bit raw files; default is 16 bit
for file = BeginPoint:(handles.PreviewNum + BeginPoint)
    
    % read files
    fid = fopen([handles.DirectoryName '\' PreviewList(file).name]);
    A = fread(fid, 'uint8=>uint8');
    fclose(fid);
    % allign bits
    E = double(A(1:2:end));
    F = double(A(2:2:end));
    G = 64*E+F/4;
    if eight_bit == 1
        intensity{file-BeginPoint+1, 1} = reshape(E, [640 480]);
    elseif eight_bit == 0
        intensity{file-BeginPoint+1, 1} = reshape(G, [640 480])';
    end
    
    % Write Files
    %     imwrite(uint16(intensity), [folder_structure '\TIFF\' current_folder '\' PreviewList(file).name '.tiff'],...
    %         'Compression', 'none');
    intensity{file-BeginPoint+1, 1} = intensity{file-BeginPoint+1, 1} - intensity0;
    
end
handles.intensity = intensity;

% Display the first one and store the graphics handle to the imshow object
handles.image = imshow(handles.intensity{1}, 'Parent', handles.axes1);

% Update the slider to accomodate all of the images
set(handles.ImageSlider, 'Min', 1, 'Max', (handles.PreviewNum + 1), ...
    'SliderStep', [1 1]/(handles.PreviewNum+1 - 1), 'Value', 1)

handles.BeginPoint = BeginPoint+handles.PreviewNum;
guidata(hObject, handles);


function ACfrequency_Callback(hObject, eventdata, handles)
% hObject    handle to ACfrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ACfrequency as text
%        str2double(get(hObject,'String')) returns contents of ACfrequency as a double

% get 'AC frequency (Hz)'
ACfrequency = str2double(get(hObject, 'String'));
handles.ACfrequency = ACfrequency;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ACfrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ACfrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function expName_Callback(hObject, eventdata, handles)
% hObject    handle to expName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expName as text
%        str2double(get(hObject,'String')) returns contents of expName as a double

% get 'Experiment No.'
expName = get(hObject, 'String');
handles.expName = expName;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function expName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ACtime_Callback(hObject, eventdata, handles)
% hObject    handle to ACtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ACtime as text
%        str2double(get(hObject,'String')) returns contents of ACtime as a double


% get 'AC time (s)'
ACtime = str2double(get(hObject, 'String'));
handles.ACtime = ACtime;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ACtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ACtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CreateROI.
function CreateROI_Callback(hObject, eventdata, handles)
% hObject    handle to CreateROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Create ROI
axes(handles.axes1);
roi = imrect;
mask = createMask(roi);
handles.mask = mask;
guidata(hObject, handles);

TimeInterval = handles.TimeInterval;
Fs = handles.Fs;
c = length(find(mask(:)~=0));

BeginPoint = 1;
FileList = dir([handles.DirectoryName '\*.raw']);
AverIntensity = zeros(length(FileList), 1);
Amp = zeros(fix(length(FileList)/Fs), 1);

if length(FileList) <= Fs*TimeInterval
    pause(TimeInterval);
end


intensity0 = handles.intensity0;
IntensityOfROI = zeros(BeginPoint+Fs*TimeInterval, 1);
while BeginPoint < length(FileList)
    EndPoint = BeginPoint+Fs*TimeInterval;
    
    if EndPoint > length(FileList)
        EndPoint = length(FileList);
    end
    
    eight_bit = 0; % select 8 bit  or 16 bit raw files; default is 16 bit
    for file = BeginPoint:EndPoint
        
        % read files
        fid = fopen([handles.DirectoryName '\' FileList(file).name]);
        A = fread(fid, 'uint8=>uint8');
        fclose(fid);
        % allign bits
        E = double(A(1:2:end));
        F = double(A(2:2:end));
        G = 64*E+F/4;
        if eight_bit == 1
            intensity = reshape(E, [640 480]);
        elseif eight_bit == 0
            intensity = reshape(G, [640 480])';
        end
        
        % Intensity of ROI (axes2)
        temp = (intensity-intensity0).*mask;
        IntensityOfROI(file, 1) = sum(temp(:))/c;
        AverIntensity(file, 1) = IntensityOfROI(file, 1);
        
        % Aver(axes4)
        
    end
    
    % Intensity of ROI (axes2)
    axes(handles.axes2);
    Lm = length(IntensityOfROI);
    handles.roiIntensity_Plot = plot((1:Lm)', IntensityOfROI);
    xlim([0 Lm])
    xlabel('Frames in interval')
    ylabel('Intensity (a.u.)')
    
    % Amplitude spectrum (axes3)
    num = fix(EndPoint/Fs);
    Y = fft(IntensityOfROI(BeginPoint:EndPoint));
    L = EndPoint - BeginPoint + 1;
    P2 = abs(Y/L);
    % Amp(num, 1) = max(2*P2(2:end-1));
    P1 = P2(1:ceil(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:ceil(L/2))/L;
    axes(handles.axes3);
    plot(f, P1)
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlim([1 fix(Fs/2)])
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    
    % Average Intensity vs time(axes4)
    axes(handles.axes4);
    LocMax = local_maximums(IntensityOfROI(BeginPoint:EndPoint));
    LocMin = local_minimums(IntensityOfROI(BeginPoint:EndPoint));
    MeanLocMax = mean(LocMax(:));
    MeanLocMin = mean(LocMin(:));
    Amp(num, 1) = abs(L*log(MeanLocMax/MeanLocMin));
    TestTime = (1:length(Amp))';
    plot(TestTime, Amp, '.');
    xlim([0 handles.ACtime])
    xlabel('t (s)')
    ylabel('Amplitude (nm)')
    
    BeginPoint = BeginPoint+Fs*TimeInterval;
    FileList = dir([handles.DirectoryName '\*.raw']);
    
    pause(0.05)
end


handles.AverIntensity = AverIntensity;
guidata(hObject, handles);
handles.Amp = Amp;
guidata(hObject, handles);
handles.TestTime = TestTime;
guidata(hObject, handles);

function op=local_maximums(s)
s1=s(1:end-2);
s2=s(2:end-1);
s3=s(3:end);
% maximums:
op=1+find((s1<=s2)&(s2>=s3));

function op=local_minimums(s)
s1=s(1:end-2);
s2=s(2:end-1);
s3=s(3:end);
% minimums:
op=1+find((s1>=s2)&(s2<=s3));

% --- Executes on button press in Delete_Save.
function Delete_Save_Callback(hObject, eventdata, handles)
% hObject    handle to Delete_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% delete ROI and save


AverIntensity = handles.AverIntensity;
Amp = handles.Amp;
TestTime = handles.TestTime;
expName = handles.expName;

[folder_structure, current_folder] = fileparts(handles.DirectoryName);
if isempty(TestTime) > 0
    mkdir([folder_structure '\MAT']);
    mkdir([folder_structure '\MAT\' current_folder  ]);
end

h = getframe(gcf);
GUIsaved = [folder_structure '\MAT\' handles.expName '_GUIsaved.tif'];
imwrite(h.cdata, GUIsaved);

SavePath = [folder_structure '\MAT\' expName '.mat'];
save(SavePath, 'AverIntensity', 'Amp', 'TestTime');


% --- Executes on button press in ZoomOn.
function ZoomOn_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.hZoom)
    handles.hZoom = zoom;
    guidata(hObject, handles);
end
if ~strcmp(get(handles.hZoom,'Enable'), 'on')
    set(handles.hZoom, 'Enable', 'on');
end


% --- Executes on button press in ZoomOff.
function ZoomOff_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.hZoom)
    set(handles.hZoom, 'Enable', 'off');
end
