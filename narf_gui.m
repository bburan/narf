function varargout = narf_gui(varargin)
% NARF_GUI MATLAB code for narf_gui.fig
%      NARF_GUI, by itself, creates a new NARF_GUI or raises the existing
%      singleton*.
%
%      H = NARF_GUI returns the handle to a new NARF_GUI or the handle to
%      the existing singleton*.
%
%      NARF_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NARF_GUI.M with the given input arguments.
%
%      NARF_GUI('Property','Value',...) creates a new NARF_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before narf_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to narf_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help narf_gui

% Last Modified by GUIDE v2.5 25-Oct-2012 10:04:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @narf_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @narf_gui_OutputFcn, ...
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


% --- Executes just before narf_gui is made visible.
function narf_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to narf_gui (see VARARGIN)

% Choose default command line output for narf_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes narf_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = narf_gui_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USEFUL FUNCTIONS FOR MORE THAN ONE WIDGET

function update_the_data_plots(hObject, eventdata, handles)
global TIME STIM RESPAVG;
% Plot the stimulus
axes(handles.stim_view_axes); cla;
nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
n = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});
plot(TIME, STIM(n,:), 'k-');
% Plot the response
axes(handles.resp_view_axes); cla;
bar(TIME, RESPAVG(n,:), 0.01,'k');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA SELECTION WIDGETS

function select_training_set_button_Callback(hObject, eventdata, handles)
% Create some globals in which to load the stimulus and response
global STIM RESP TIME RESPAVG;

% Use David's Nifty DB file chooser to choose the file data
fd = dbchooserawfile(0,'Choose file to sort');
if isempty(fd)
    STIM = [];
    RESP = [];
    set(handles.training_set, 'String', 'none');
    guidata(hObject, handles);               % Write the change
    drawnow;                                 % Flush the event queue
    return
end

% Q: WHY IS IT SO HARD TO USE CELLDB TO GET A FILE!?!?
% A: Because it really limits what you can query by!
[cfd, cellids, cellfileids] = dbgetscellfile('rawid', fd.rawid);

% TODO: Is this right?
% Unfortunately, we must hack around the fact that there is typically
% duplicate cellfiledescriptors (cfds) returned by the above, but the
% one we want (I think) is the one with respfiletype=1, because
% (I think) that means it has been spike-sorted already?
% BUT I'VE BEEN WRONG BEFORE!
index = find([cfd.respfiletype] == 1);

% If there is more than one cell file returned, complain
if size(index) ~= 1
    error('BAPHY gave me multiple files and I cannot find the one I need');
end

options.includeprestim = 1;
options.unit     = cfd(index).unit;
options.channel  = fd.channel;
options.rasterfs = 100000; 

fprintf('Loading stimulus file: %s%s\n', cfd(index).stimpath, cfd(index).stimfile);
stimfile = [cfd(index).stimpath cfd(index).stimfile];
stim     = loadstimfrombaphy(stimfile, [], [], 'wav', options.rasterfs, 1, 0, options.includeprestim);

fprintf('Loading response file: %s/%s\n', cfd(index).path, cfd(index).respfile);
respfile = [cfd(index).path cfd(index).respfile];
[resp, tags] = loadspikeraster(respfile, options);

[d1 d2 d3] = size(stim);
if d1 ~= 1
    disp([d1 d2 d3]);
    error('Stimulus matrix must initially have size 1xLxN, where L=length in samples, N=num of stimuli');
end

STIM = reshape(stim, [d3 d2]);  % Remove the irrelevant first dimension
RESP = permute(resp, [3 1 2]);  % Rearrange to match
rsp = sum(RESP,3); 
[e1 e2 e3] = size(rsp);
RESPAVG = reshape(rsp, [e1 e2]); % Averaged response across all trials
TIME = (1/options.rasterfs).*[1:d2]';

% TODO: Check that STIM, RESP, RESPAVG, TIME are all the proper sizes

% Update other parts of the GUI, starting with the text box
set(handles.training_set, 'String', cellids(index));   

% Update the popup for selecting a particular stimuli to graph
c = {};
for i=1:e1;
    c{i} = sprintf('%d',i);
end
set(handles.selected_stimuli_popup, 'String', char(c)); 
set(handles.selected_stimuli_popup, 'Value', 1);

update_the_data_plots();

% DELETE ME?
%guidata(hObject, handles);               % Write the change
%drawnow;                                 % Flush the event queue

function refresh_stim_view_button_Callback(hObject, eventdata, handles)
update_the_data_plots(hObject, eventdata, handles);

function select_test_set_button_Callback(hObject, eventdata, handles)

function training_set_info_listbox_CreateFcn(hObject, eventdata, handles)

function training_set_info_listbox_Callback(hObject, eventdata, handles)

function selected_stimuli_popup_CreateFcn(hObject, eventdata, handles)

function selected_stimuli_popup_Callback(hObject, eventdata, handles)
update_the_data_plots(hObject, eventdata, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPROCESSING WIDGETS

function apply_prefilter_button_Callback(hObject, eventdata, handles)

% Get the low-pass and high-pass band limits
lp = eval(get(findobj(gcf,'Tag','low_freqs'), 'String'));
hp = eval(get(findobj(gcf,'Tag','high_freqs'), 'String'));
sf = eval(get(findobj(gcf,'Tag','smoothing_frqs'), 'String'));
    
%ah = findobj(gcf,'Tag','stim_view_axes');      % Update stimulus view
%hold on;
%for i = 1:length()
%    plot(ah, TIME, filter(ellip(4,0.5,20,[lp(i) hp(i)]),STIM), 'k-');
%end
%hold off;
%set(ah,'Tag','stim_view_axes');

function popupmenu1_CreateFcn(hObject, eventdata, handles)

function popupmenu1_Callback(hObject, eventdata, handles)

function preproc_filter_popup_CreateFcn(hObject, eventdata, handles)
% TODO: Get a list of all the preprocessing filters that could work
% hObject.String = char({'Option1', 'Option2', 'Option3'});

function preproc_filter_popup_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
disp(get(hObject,'String'));
fn_to_run = contents{get(hObject,'Value')};
% 
% function edit_pre_params_button_Callback(hObject, eventdata, handles)
% % TODO: Rewrite this with axes(handles.axes1);cla;
% % popup_sel_index = get(handles.popupmenu1, 'Value');
% disp('pushbutton 1 was pressed');
% ah = findobj(gcf,'Tag','plotarea_preprocessed');
% x=rand(randi(10+20,1),4);
% plot(ah, x);
% set(ah,'Tag','plotarea_preprocessed'); % Don't know why it needs this.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CLASS WIDGETS

function fir_history_CreateFcn(hObject, eventdata, handles)

function fir_history_Callback(hObject, eventdata, handles)
str2double(get(hObject,'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STOCHASTICITY WIDGETS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMIZATION WIDGETS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORMANCE WIDGETS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fir_filter_popup_CreateFcn(hObject, eventdata, handles)

function fir_filter_popup_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns fir_filter_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fir_filter_popup


function smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothing as text
%        str2double(get(hObject,'String')) returns contents of smoothing as a double


% --- Executes during object creation, after setting all properties.
function smoothing_CreateFcn(hObject, eventdata, handles)


function edit8_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dimensions_Callback(hObject, eventdata, handles)
% hObject    handle to dimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dimensions as text
%        str2double(get(hObject,'String')) returns contents of dimensions as a double


% --- Executes during object creation, after setting all properties.
function dimensions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function low_freqs_Callback(hObject, eventdata, handles)
% hObject    handle to low_freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_freqs as text
%        str2double(get(hObject,'String')) returns contents of low_freqs as a double


% --- Executes during object creation, after setting all properties.
function low_freqs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function high_freqs_Callback(hObject, eventdata, handles)
% hObject    handle to high_freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of high_freqs as text
%        str2double(get(hObject,'String')) returns contents of high_freqs as a double


% --- Executes during object creation, after setting all properties.
function high_freqs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to high_freqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in model_view_popup.
function model_view_popup_Callback(hObject, eventdata, handles)
% hObject    handle to model_view_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_view_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_view_popup


% --- Executes during object creation, after setting all properties.
function model_view_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_view_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_model_button.
function refresh_model_button_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_model_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in data_view_popup.
function data_view_popup_Callback(hObject, eventdata, handles)
% hObject    handle to data_view_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_view_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_view_popup


% --- Executes during object creation, after setting all properties.
function data_view_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_view_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_data_button.
function refresh_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in preprocessing_view_popup.
function preprocessing_view_popup_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessing_view_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns preprocessing_view_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from preprocessing_view_popup


% --- Executes during object creation, after setting all properties.
function preprocessing_view_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to preprocessing_view_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_preproc_button.
function refresh_preproc_button_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_preproc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function smoothing_frqs_Callback(hObject, eventdata, handles)
% hObject    handle to smoothing_frqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothing_frqs as text
%        str2double(get(hObject,'String')) returns contents of smoothing_frqs as a double


% --- Executes during object creation, after setting all properties.
function smoothing_frqs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothing_frqs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
