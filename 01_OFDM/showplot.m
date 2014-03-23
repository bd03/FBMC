function varargout = showplot(varargin)
% SHOWPLOT MATLAB code for showplot.fig
%      SHOWPLOT, by itself, creates a new SHOWPLOT or raises the existing
%      singleton*.
%
%      H = SHOWPLOT returns the handle to a new SHOWPLOT or the handle to
%      the existing singleton*.
%
%      SHOWPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWPLOT.M with the given input arguments.
%
%      SHOWPLOT('Property','Value',...) creates a new SHOWPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before showplot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to showplot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help showplot

% Last Modified by GUIDE v2.5 21-Mar-2014 02:12:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showplot_OpeningFcn, ...
                   'gui_OutputFcn',  @showplot_OutputFcn, ...
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

% --- Executes just before showplot is made visible.
function showplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showplot (see VARARGIN)

% Choose default command line output for showplot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%initial plot
pushbutton2_Callback(hObject, eventdata, handles)
pushbutton4_Callback(hObject, eventdata, handles)
update_plot(handles)


% UIWAIT makes showplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = showplot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1
update_plot(handles)

% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2
update_plot(handles)

% --- Executes on button press in togglebutton3.
function togglebutton3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton3
update_plot(handles)

% --- Executes on button press in togglebutton4.
function togglebutton4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton4
update_plot(handles)

% --- Executes on button press in togglebutton20.
function togglebutton20_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton20
update_plot(handles)

% --- Executes on button press in togglebutton21.
function togglebutton21_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton21
update_plot(handles)

% --- Executes on button press in togglebutton22.
function togglebutton22_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton22
update_plot(handles)

% --- Executes on button press in togglebutton23.
function togglebutton23_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton23
update_plot(handles)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_array=[handles.togglebutton20,handles.togglebutton21, ...
    handles.togglebutton22, handles.togglebutton23,handles.togglebutton24];
for i=1:length(handles_array)
    set(handles_array(i),'Value',0)
end
update_plot(handles);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_array=[handles.togglebutton20,handles.togglebutton21, ...
    handles.togglebutton22, handles.togglebutton23,handles.togglebutton24];
for i=1:length(handles_array)
    set(handles_array(i),'Value',1)
end

update_plot(handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_array =[handles.togglebutton1,handles.togglebutton13,...
    handles.togglebutton14,handles.togglebutton15,handles.togglebutton16,...
    handles.togglebutton17, handles.togglebutton18,handles.togglebutton19];
for i=1:length(handles_array)
    set(handles_array(i),'Value',0)
end
update_plot(handles);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_array =[handles.togglebutton1,handles.togglebutton13,...
    handles.togglebutton14,handles.togglebutton15,handles.togglebutton16,...
    handles.togglebutton17, handles.togglebutton18,handles.togglebutton19];
for i=1:length(handles_array)
    set(handles_array(i),'Value',1)
end
update_plot(handles);

% --- Executes on button press in togglebutton13.
function togglebutton13_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton13
update_plot(handles)

% --- Executes on button press in togglebutton14.
function togglebutton14_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton14
update_plot(handles)

% --- Executes on button press in togglebutton15.
function togglebutton15_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton15
update_plot(handles)

% --- Executes on button press in togglebutton16.
function togglebutton16_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton16
update_plot(handles)

% --- Executes on button press in togglebutton17.
function togglebutton17_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton17
update_plot(handles)

% --- Executes on button press in togglebutton18.
function togglebutton18_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton18
update_plot(handles)

% --- Executes on button press in togglebutton19.
function togglebutton19_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton19
update_plot(handles)

% --- Executes on button press in togglebutton24.
function togglebutton24_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton24
update_plot(handles)

function update_plot(handles)
% toggle1: M=4
% toggle13-19: M=8,16,32,64,128,256,512
% toggle20-24: 4-16-64-128-256-QAM
% toggle button handles belonging to different M values
handles_M =[handles.togglebutton1,handles.togglebutton13,...
    handles.togglebutton14,handles.togglebutton15,handles.togglebutton16,...
    handles.togglebutton17, handles.togglebutton18,handles.togglebutton19];

% toggle button handles belonging to different modulation values
handles_mod=[handles.togglebutton20,handles.togglebutton21, ...
    handles.togglebutton22, handles.togglebutton23, handles.togglebutton24];

%SNR array, when integrated, this can get values directly from
%FBMC_simulation file
% load('SNR_array.mat');
% load('qam_sizes.mat');
% load('M_array.mat');
% file = ls('BER2014*');
% load(file(size(file,1),:));
try
    load('BER_archive/BER2014-3-22-20-15-Final.mat');
catch
    error('BER data could not be found.')
end

try
    load('BER_archive/CONF2014-3-22-20-15-Final.mat');
catch
    error('CONF data could not be found.')
end

% process conf
M_array = conf(1).M_val;
qam_sizes = conf(1).mod_sch;
SNR_array = conf(1).SNR_val;
all_mod =[4 16 64 128 256];
all_M =[4 8 16 32 64 128 256 512];

%active buttons
for i=all_M
    if ~ismember(i, M_array)
        set(handles_M(log2(i)-1),'Value',0,'Enable','off');
    end
end

pp=1;
for i=all_mod
    if ~ismember(i, qam_sizes)
        set(handles_mod(find(all_mod==i)),'Value',0,'Enable','off');
    end
end

% retrieve data
data=zeros(length(M_array)*length(qam_sizes),length(SNR_array));
for qw=1:length(M_array)
    for zx=1:length(qam_sizes);
        data(length(qam_sizes)*(qw-1)+zx,:)=BER(qw,zx,SNR_array+1);
    end
end

% plot data 
% data contains BER for (16) different SNR values and oriented in following 
% fashion:
% M=4/4-QAM/SNR=0 M=4/4-QAM/SNR=1 M=4/4-QAM/SNR=2 ..........
% M=4/16-QAM/SNR=0    ......    ......    ......    ......    ......    
% ......    ......    ......    ......    ......    ......    ......    
% ......    ......    ......    ......    ......    ......    ......    1
% M=512/128-QAM/SNR=0............................M=512/128-QAM/SNR=15

% line modifiers
markers= ['+' 'o' '*' 's' '^'];
lines = ['-' ':' '-.' '--'];
colors = ['r' 'g' 'b' 'm' 'k'];
combined = ['--g' '-m' '--b' '-r' '--k' '-g' '--m' '-b' '--r' '-k'];

%colors and styles
flag_at_least_one = false;
for qw=M_array
    for zx=qam_sizes
        if get(handles_M(log2(qw)-1),'Value')==1 && get(handles_mod(find(all_mod==zx)),'Value')==1 && ...
                strcmp(get(handles_M(log2(qw)-1),'Enable'),'on') && strcmp(get(handles_mod(find(all_mod==zx)),'Enable'),'on') 
            %disp(sprintf('qw%d zx%d',qw,zx))
            modifier = strcat(lines(mod(log2(qw)-1-1,2)+1),colors(mod(log2(qw)-1,5)+1),markers(mod(find(all_mod==zx)-1,5)+1));
            semilogy(SNR_array,data(length(qam_sizes)*(find(M_array==qw)-1)+find(qam_sizes==zx),:),modifier,...
                'LineWidth',1.5,'MarkerSize',8);
            flag_at_least_one = true;
            hold on
            grid on
        end
    end
end
xlabel('SNR (dB)');
ylabel('BER');
hold off

if ~flag_at_least_one
    semilogy(1,1);
end
