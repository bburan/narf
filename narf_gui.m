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

% Last Modified by GUIDE v2.5 26-Oct-2012 14:59:02

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

function c = pickcolor(n)
m = mod(n,7);
switch m
    case 0 
        c='k-';
    case 1 
        c='b-';
    case 2 
        c='g-';
    case 3 
        c='r-';
    case 4 
        c='c-';
    case 5
        c='m-';
    case 6
        c='y-';
end

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

function update_prefilter_plots(hObject, eventdata, handles)
global TIME STIM PREFILTERS RESPAVG PREFILTEREDSTIM;

nvals = cellstr(get(handles.preprocessing_view_popup, 'String'));
plottype = nvals{get(handles.preprocessing_view_popup, 'Value')};

nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
n = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});

axes(handles.preproc_view_axes); cla;
hold on;
switch plottype
    case 'Freq. Response'
        for i=length(PREFILTERS)
            [sos,g] = tf2sos(PREFILTERS{i}{1}, PREFILTERS{i}{2}); 
            Hd = dfilt.df2tsos(sos,g);  % Create a dfilt object
            h = fvtool(Hd);              % Plot magnitude response
            set(h,'Analysis','magnitude') % Freq Magnitude response
        end
    case 'Filtered Stimulus'
        for i=length(PREFILTERS)  
            ps = cell2mat(PREFILTEREDSTIM{i});
            plot(TIME, ps(n,:), pickcolor(i));
        end
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA SELECTION WIDGETS

function select_training_set_button_Callback(hObject, eventdata, handles)

% Create some globals in which to load the stimulus and response
global STIM RESP TIME RESPAVG SAMPFREQ;

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
SAMPFREQ = options.rasterfs; 

fprintf('Loading stimulus file: %s%s\n', cfd(index).stimpath, cfd(index).stimfile);
stimfile = [cfd(index).stimpath cfd(index).stimfile];
stim     = loadstimfrombaphy(stimfile, [], [], 'wav', options.rasterfs, 1, 1, options.includeprestim);

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

update_the_data_plots(hObject, eventdata, handles);

function refresh_stim_view_button_Callback(hObject, eventdata, handles)
update_the_data_plots(hObject, eventdata, handles);

function select_test_set_button_Callback(hObject, eventdata, handles)

function training_set_info_listbox_CreateFcn(hObject, eventdata, handles)
function training_set_info_listbox_Callback(hObject, eventdata, handles)

function selected_stimuli_popup_CreateFcn(hObject, eventdata, handles)
function selected_stimuli_popup_Callback(hObject, eventdata, handles)
update_the_data_plots(hObject, eventdata, handles);

function data_view_popup_Callback(hObject, eventdata, handles)
function data_view_popup_CreateFcn(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPROCESSING WIDGETS

function preproc_filter_popup_CreateFcn(hObject, eventdata, handles)
% TODO: On startup, find all the preprocessing filters in a directory

function preproc_filter_popup_Callback(hObject, eventdata, handles)
%contents = cellstr(get(hObject,'String'));
%disp(get(hObject,'String'));
%fn_to_run = contents{get(hObject,'Value')};

function low_freqs_CreateFcn(hObject, eventdata, handles)
function low_freqs_Callback(hObject, eventdata, handles)
function high_freqs_CreateFcn(hObject, eventdata, handles)
function high_freqs_Callback(hObject, eventdata, handles)
function smoothing_frqs_Callback(hObject, eventdata, handles)
function smoothing_frqs_CreateFcn(hObject, eventdata, handles)

function preprocessing_view_popup_CreateFcn(hObject, eventdata, handles)
function preprocessing_view_popup_Callback(hObject, eventdata, handles)
update_prefilter_plots(hObject, eventdata, handles);

function apply_prefilter_button_Callback(hObject, eventdata, handles)
global STIM SAMPFREQ PREFILTERS PREFILTEREDSTIM;
[d1, d2] = size(STIM);
% Get the low-pass and high-pass band limits
lp = eval(get(handles.low_freqs, 'String'));
hp = eval(get(handles.high_freqs, 'String'));
sf = eval(get(handles.smoothing_frqs, 'String'));
% Check that the inputs are valid
if ~(isvector(lp) & isvector(hp) & isvector(sf))
    set(handles.prefilter_status, 'String', 'ERROR: ALL PREFILTER SETTINGS MUST BE VECTORS AND OF THE SAME LENGTH!');
    return;
end
n = length(lp); % How many filters are there?
PREFILTERS={};  % Make the bank of filters
for i = 1:n
    [B,A] = ellip(4,0.5,50,[lp(i)/SAMPFREQ*2, hp(i)/SAMPFREQ*2]);   
    PREFILTERS{i} = {B,A};
end
% Filter the data
set(handles.prefilter_status, 'String', 'Filtering...');
PREFILTEREDSTIM={};
for i = 1:n
    pfstmp1 = filter(PREFILTERS{i}{1}, PREFILTERS{i}{2}, STIM);
    pfstmp2 = abs(pfstmp1);
    % OPTIONAL LOW_PASS SMOOTHING ALGORITHM
    %     [B,A] = ellip(6,3,50,[sf(i)/SAMPFREQ*2]);
    %     pfs = filter(B,A, pfstmp2);
    PREFILTEREDSTIM{i} = {pfstmp2};
end
% Plot either the data or the frequency response
set(handles.prefilter_status, 'String', 'Plotting...');
update_prefilter_plots(hObject, eventdata, handles);
set(handles.prefilter_status, 'String', 'Done.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CLASS WIDGETS

function model_class_popup_CreateFcn(hObject, eventdata, handles)
function model_class_popup_Callback(hObject, eventdata, handles)
% TODO: In the future, create 

function fir_history_CreateFcn(hObject, eventdata, handles)
function fir_history_Callback(hObject, eventdata, handles)
% TODO: Initialize the FIR filter_dimensions
str2double(get(hObject,'String'));

function model_view_popup_CreateFcn(hObject, eventdata, handles)
function model_view_popup_Callback(hObject, eventdata, handles)

function plot_model_button_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMIZATION WIDGETS

function sampling_algorithm_popup_CreateFcn(hObject, eventdata, handles)
function sampling_algorithm_popup_Callback(hObject, eventdata, handles)

function performance_metric_popup_CreateFcn(hObject, eventdata, handles)
function performance_metric_popup_Callback(hObject, eventdata, handles)

function stochasticity_popup_CreateFcn(hObject, eventdata, handles)
function stochasticity_popup_Callback(hObject, eventdata, handles)

function termination_condition_popup_CreateFcn(hObject, eventdata, handles)
function termination_condition_popup_Callback(hObject, eventdata, handles)

function termination_iterations_CreateFcn(hObject, eventdata, handles)
function termination_iterations_Callback(hObject, eventdata, handles)

function fit_model_button_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORMANCE WIDGETS

% Nothing so far!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
