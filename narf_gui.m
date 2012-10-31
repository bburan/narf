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

% ADD some things to the path!
NARFHOME = '/home/ivar/matlab/narf'
addpath([NARFHOME filesep 'utils'], ...
        [NARFHOME filesep 'samplers']);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA SELECTION WIDGETS

function update_the_data_plots(hObject, eventdata, handles)
global TIME STIM RESPAVG SAMPFREQ;

% Plot the stimulus
axes(handles.stim_view_axes); cla;

nsel = cellstr(get(handles.data_view_popup, 'String'));
plottype = nsel{get(handles.data_view_popup, 'Value')};
nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
n = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});
switch plottype
    case 'Time Series View'
        plot(TIME, STIM(n,:), 'k-');
    case 'Spectrogram View'
        % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
        nwin = 4048;
        logfsgram(STIM(n,:)', nwin, SAMPFREQ, [], [], 500, 12);
        caxis([-20,40]);
        %colorbar;
end

% Plot the response
axes(handles.resp_view_axes); cla;
bar(TIME, RESPAVG(n,:), 0.01,'k');

%------------------------------------------------------------------------
function select_training_set_button_Callback(hObject, eventdata, handles)
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

[cfd, cellids, cellfileids] = dbgetscellfile('rawid', fd.rawid);

% TODO: Replace this file with the right channel and unit?
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

STIM = permute(stim, [3 2 1]);  % Remove the irrelevant first dimension
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
%------------------------------------------------------------------------

function refresh_stim_view_button_Callback(hObject, eventdata, handles)
update_the_data_plots(hObject, eventdata, handles);

function select_test_set_button_Callback(hObject, eventdata, handles)

function training_set_info_listbox_CreateFcn(hObject, eventdata, handles)
function training_set_info_listbox_Callback(hObject, eventdata, handles)

function selected_stimuli_popup_CreateFcn(hObject, eventdata, handles)
function selected_stimuli_popup_Callback(hObject, eventdata, handles)
update_the_data_plots(hObject, eventdata, handles);

function data_view_popup_CreateFcn(hObject, eventdata, handles)
function data_view_popup_Callback(hObject, eventdata, handles)
update_the_data_plots(hObject, eventdata, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPROCESSING WIDGETS

function update_prefilter_plots(handles)
global TIME STIM PREFILTERS RESPAVG PREFILTEREDSTIM;

nvals = cellstr(get(handles.preprocessing_view_popup, 'String'));
plottype = nvals{get(handles.preprocessing_view_popup, 'Value')};

nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
n = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});
    
axes(handles.preproc_view_axes); cla;
hold on;
switch plottype
    case 'Freq. Response'
        for i=1:length(PREFILTERS)
            [sos,g] = tf2sos(PREFILTERS{i}{1}, PREFILTERS{i}{2}); 
            Hd = dfilt.df2tsos(sos,g);  % Create a dfilt object
            h = fvtool(Hd);              % Plot magnitude response
            set(h,'Analysis','magnitude') % Freq Magnitude response
        end
    case 'Filtered Stimulus'
        for i=1:length(PREFILTERS)  
            ps = cell2mat(PREFILTEREDSTIM{i});
            plot(TIME, ps(i,:), pickcolor(i));
        end
end
hold off;

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
update_prefilter_plots(handles);


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
    pfstmp1 = filter(PREFILTERS{i}{1}, PREFILTERS{i}{2}, STIM,[],2);
    pfstmp2 = abs(pfstmp1);
    % OPTIONAL LOW_PASS SMOOTHING ALGORITHM
    if sf ~= 0
        [B,A] = ellip(6,3,50,[sf(i)/SAMPFREQ*2]);
        pfs = filter(B,A, pfstmp2,[],2);
    else
        pfs = pfstmp2;
    end
    PREFILTEREDSTIM{i} = {pfs};
end

% Plot either the data or the frequency response
set(handles.prefilter_status, 'String', 'Plotting...');
update_prefilter_plots(handles);
set(handles.prefilter_status, 'String', 'Done.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CLASS WIDGETS
function update_fir_model(handles)
global PREFILTERS FIRHIST FIRBINSIZE FIRCOEFS;
FIRHIST = str2num(get(handles.fir_history,'String'));
FIRBINSIZE = str2num(get(handles.bin_size,'String'));
d = length(PREFILTERS);
l = ceil(FIRHIST/FIRBINSIZE);
FIRCOEFS = zeros(d,l);

function downsample_stimresp(stim, resp, sampfrq, binsize)
% Downsample (DS_) the stimulus and response 
global PREFILTEREDSTIM RESPAVG SAMPFREQ DS_STIM DS_RESP;

c = length(PREFILTEREDSTIM); % Number of filter inputs
[d,l] = size(RESPAVG); % d=num of tones, l=num of samples

% TODO: The following algorithm may miss a few elements in the last bin
%    It should be fixed so that they are included.
 nagg = ceil(SAMPFREQ / (1000*FIRBINSIZE)); % Number of samples to aggregate
 nl = floor(l/nagg); % The number of samples in the new system
 BINNEDRESP = zeros(d, nl);
 BINNEDSTIM = zeros(d, nl, c);
 
 disp('binning stim and resp');
% Bin the stimulus and response signals
for i = 1:d
    for j = 1:nl
        DS_RESP(i,j) = mean(RESPAVG(i,(nagg*(j-1)+1):nagg*j));
        for k = 1:c
            tmp = cell2mat(PREFILTEREDSTIM{k});
            DS_STIM(i,j,k) = mean(tmp(i,(nagg*(j-1)+1):nagg*j));
        end
    end
end


function plot_model(handles)
global FIRCOEFS DSTIME ;
nvals = cellstr(get(handles.model_view_popup, 'String'));
plottype = nvals{get(handles.model_view_popup, 'Value')};
[d, l] = size(FIRCOEFS);

axes(handles.model_view_axes); cla;
hold on;
switch plottype
    case 'FIR Shape'
        for i=1:d
            stem([1:l], FIRCOEFS(i,:), pickcolor(i));
        end
    case 'Stim/Resp/Pred'
        plot(
        for i=1:d
            % plot the down sampled stim, response, and model prediction
            plot(TIME, ones(1,length(TIME)), pickcolor(i));
            
        end
end
hold off;

function model_class_popup_CreateFcn(hObject, eventdata, handles)
function model_class_popup_Callback(hObject, eventdata, handles)
% TODO: In the future choose the model class here 

function fir_history_CreateFcn(hObject, eventdata, handles)
function fir_history_Callback(hObject, eventdata, handles)
update_fir_model(handles);

function model_view_popup_CreateFcn(hObject, eventdata, handles)
function model_view_popup_Callback(hObject, eventdata, handles)
plot_model(handles);

function plot_model_button_Callback(hObject, eventdata, handles)
plot_model(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMIZATION WIDGETS

function x = pak_FIRCOEFS(FC)
x = FC(:);

function FC = unpak_FIRCOEFS(x)
global PREFILTERS FIRHIST FIRBINSIZE;
d = length(PREFILTERS);
l = ceil(FIRHIST/FIRBINSIZE);
FC = reshape(x, d, l);

function z = binned_correlation(x)
% Returns the correlation of the binned stimulus and binned response
global RESPAVG PREFILTEREDSTIM FIRBINSIZE SAMPFREQ FIRCOEFS;

%fprintf('binned_correlation(x):');
%disp(x');

unpak_FIRCOEFS(x);

 
% Apply the FIR filter to the Binned Stimulus to get the model prediction
% Since it is linear, the prediction is just the sum of the filters
% We assume that there are no second order terms combining elements of both filters
disp('Applying FIR filter');
[n,l] = size(FIRCOEFS);
for i = 1:n
    if i == 1
            
        
        
        
        BINNEDPRED = filter(FIRCOEFS(i,:), , BINNEDSTIM, [],2);
    else
        BINNEDPRED = BINNEDPRED + filter(FIRCOEFS(1,:), FIRCOEFS(2,:), BINNEDSTIM, [],2);
    end
end

% Flatten the prediction and response and compare using correlation
R = corrcoef(BINNEDPRED, BINNEDRESP);
z = R(2,1);

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
% get the number of iterations
niters = str2num(get(handles.termination_iterations, 'String'));

% Get the starting score of the current filter
x_0 = pak_FIRCOEFS();

stepsize = 1.0;
[x_bst, s_bst] = boosting(x_0, @binned_correlation, @(n,x,s)(n > niters), stepsize);
unpak_FIRCOEFS(x_bst);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REPORTING WIDGETS

% Nothing so far!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
