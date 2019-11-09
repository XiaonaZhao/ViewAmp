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

% Last Modified by GUIDE v2.5 08-Nov-2019 21:58:11

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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[TifFile, TifPath] = uigetfile('*.tiff', 'Multiselect', 'off', 'Load a tif File');
OneShot = imread(fullfile(TifPath, TifFile));
handles.OneShot = OneShot;
guidata(hObject, handles);
axes(handles.axes1);
imshow(OneShot, 'DisplayRange', [], 'InitialMagnification', 'fit');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
roi = imrect;
mask = createMask(roi);
handles.mask = mask;
guidata(hObject, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% delete ROI
OneShot = handles.OneShot;
imshow(OneShot, 'DisplayRange', [], 'InitialMagnification', 'fit', handles.axes1);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TimeInterval = handles.TimeInterval;
Fs = handles.Fs;
mask = handles.mask;
c = length(find(mask(:)~=0));

DirectoryName = uigetdir();
handles.DirectoryName = DirectoryName;
guidata(hObject,handles);

BeginPoint = 1;
FileList = dir([DirectoryName '\*.raw']);
AverIntensity = zeros(length(FileList), 1);
Amp = zeros(fix(length(FileList)/Fs), 1);

if length(FileList) <= Fs*TimeInterval
    pause(TimeInterval);
end

while BeginPoint < length(FileList)
    EndPoint = BeginPoint+Fs*TimeInterval;
    
    if EndPoint > length(FileList)
        EndPoint = length(FileList);
    end
    IntensityOfROI = zeros(EndPoint -BeginPoint+1, 1);
    
    eight_bit = 0; % select 8 bit  or 16 bit raw files; default is 16 bit
    for file = BeginPoint:EndPoint
        
        % read files
        fid = fopen([DirectoryName '\' FileList(file).name]);
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
        temp = intensity.*mask;
        IntensityOfROI(file-BeginPoint+1, 1) = sum(temp(:))/c;
        AverIntensity(file, 1) = IntensityOfROI(file-BeginPoint+1, 1);
        
        % Aver(axes4)
        
    end
    
    % Intensity of ROI (axes2)
    axes(handles.axes2);
    L = length(IntensityOfROI);
    plot((1:L)', IntensityOfROI);
    xlim([0 Fs*TimeInterval])
    xlabel('Frames in interval')
    ylabel('Intensity (a.u.)')
    
    % Amplitude spectrum (axes3)
    num = fix(EndPoint/Fs);
    Y = fft(IntensityOfROI);
    P2 = abs(Y/L);
    Amp(num, 1) = max(2*P2(2:end-1));
    P1 = P2(1:ceil(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:ceil(L/2))/L;
    axes(handles.axes3);
    plot(f, P1)
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlim([1 fix(Fs/2)])
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    %     pause(0.3);
    
    % Average Intensity vs time(axes4)
    axes(handles.axes4);
    t = (1:length(Amp))';
    plot(t, Amp, '.');
    xlim([0 50])
    xlabel('t (s)')
    ylabel('Intensity Amp per sec')
    %     pause(0.3);
    
    BeginPoint = BeginPoint+Fs*TimeInterval;
    FileList = dir([DirectoryName '\*.raw']);
end

handles.AverIntensity = AverIntensity;
guidata(hObject, handles);
handles.Amp = Amp;
guidata(hObject, handles);
handles.t = t;
guidata(hObject, handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% STOP and Save plots
DirectoryName = handles.DirectoryName;
AverIntensity = handles.AverIntensity;
Amp = handles.Amp;
t = handles.t;
expName = handles.expName;

SavePath = [DirectoryName '\' expName '.mat'];
save(SavePath, 'AverIntensity', 'Amp', 't');


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% get 'Sample rate (fps)'
Fs = str2double(get(hObject, 'String'));
handles.Fs = Fs;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

% get 'AC frequency (Hz)'
ACfrequency = str2double(get(hObject, 'String'));
handles.ACfrequency = ACfrequency;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

% get 'Time interval (s)'
TimeInterval = str2double(get(hObject, 'String'));
handles.TimeInterval = TimeInterval;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double

% get 'Experiment No.'
expName = get(hObject, 'String');
handles.expName = expName;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
