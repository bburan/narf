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

% Last Modified by GUIDE v2.5 07-Nov-2012 11:38:41

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION CODE EXECUTES BEFORE GUI BECOMES VISIBLE

function narf_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to narf_gui (see VARARGIN)

% Choose default command line output for narf_gui
handles.output = hObject;
%handles.output = handles.status_textbox; % TODO: Enable this

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes narf_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Add necessary directories to NARF's path
NARFHOME = '/home/ivar/matlab/narf'
addpath([NARFHOME filesep 'utils'], ...
        [NARFHOME filesep 'samplers']);

% Create a scrollable edit box in the status panel
hEdit = uicontrol(handles.uipanel6, 'Style','edit', 'FontSize',9, ...
    'Min',0, 'Max',2, 'HorizontalAlignment','left', ...
    'Units','normalized', 'Position',[0 0 1 1], ...
    'String','GUI Initialized');

global LOG_HANDLE LOG_LENGTH LOG_BUFFER LOG_LEVEL;
LOG_LEVEL = 0;      % 0=debug, 1=informative, 2=normal, 3=warnings, 4=errors % TODO: Make top-level constant
LOG_HANDLE = hEdit;
LOG_LENGTH = 6;      % TODO: Make a top-level constant
LOG_BUFFER = cell(1,LOG_LENGTH);
for i = 1:LOG_LENGTH
    LOG_BUFFER{i} = '';
end

% Invalidate all data tables
set(handles.data_selection_table, 'Data', {});
set(handles.sampling_data_table, 'Data', {});
set(handles.perf_metric_data_table, 'Data', {});
set(handles.term_cond_data_table, 'Data', {});
set(handles.preproc_data_table, 'Data', {});
set(handles.model_data_table, 'Data', {});
set(handles.stochasticity_data_table, 'Data', {});
drawnow;

% --- Outputs from this function are returned to the command line.
function varargout = narf_gui_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USEFUL FUNCTIONS FOR MORE THAN ONE WIDGET

function printlog(varargin)
global LOG_HANDLE LOG_LENGTH LOG_BUFFER;
for i = 1:LOG_LENGTH-1                                % Shift buffer up
    LOG_BUFFER{i} = LOG_BUFFER{i+1};
end
LOG_BUFFER{LOG_LENGTH} = feval(@sprintf, varargin{:});% Add the newest data
set(LOG_HANDLE, 'String', char(LOG_BUFFER));          % Refresh the display
drawnow;

function log_dbg(varargin) % For debugging-level messages only
global LOG_LEVEL; if(LOG_LEVEL < 1)  feval(@printlog, varargin{:}); end

function log_inf(varargin) % For informative messages
global LOG_LEVEL; if(LOG_LEVEL < 2)  feval(@printlog, varargin{:}); end

function log_msg(varargin) % For regular priority messages
global LOG_LEVEL; if(LOG_LEVEL < 3)  feval(@printlog, varargin{:}); end

function log_wrn(varargin) % For warning messages
global LOG_LEVEL; if(LOG_LEVEL < 4)  feval(@printlog, varargin{:}); end

function log_err(varargin) % For error messages
global LOG_LEVEL; if(LOG_LEVEL < 5)  feval(@printlog, varargin{:}); end
error(feval(@sprintf, varargin{:}));

% ------------------------------------------------------------------------
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
% Core functionality of the system should be independent of the GUI

function run_narf(cellid, dataset_sel_fn, opt_fn, perf_metric, term_cond)
cfd = query_db(cellid);
[train_set, test_set] = select_train_test_sets(cfd); 
% Load the files
% -> Display the raw data
% Load the model and its parameters
% Initialize the model
% Prefilter the data
% -> Display the filtered data
% Downsample the filtered data
% -> Display the downsampled data
% Run the optimization on the downsampled data
% -> Display optimization performance
% -> Display parameters in GUI
% Return the best fitting parameters

% ------------------------------------------------------------------------
function cfd = query_db(cellid)
log_dbg('query_db(''%s'');', cellid);

% Returns a list of raw files with this cell ID
[cfd, cellids, cellfileids] = dbgetscellfile('cellid', cellid);

% If there is not exactly one cell file returned, throw an error.
if ~isequal(length(cellids), 1)
    log_err('BAPHY gave me zero/multiple cellids yet I only need one.');
end

% ------------------------------------------------------------------------
function [training_set, test_set] = select_train_test_sets(cfd)
log_dbg('Selecting training and test sets...');
len = length(cfd);
training_set = {};
test_set = {};
test_set_reps = 0;
parms = cell(1, len);
perfs = cell(1, len);

% Load parms, perfs. Select set with the most repetitions as test set.
for i = 1:len;
    [parms{i}, perfs{i}] = dbReadData(cfd(i).rawid);
    if (isfield(parms{i}, 'Ref_Subsets') & ...
            parms{i}.Ref_Subsets > test_set_reps)
       test_set_reps = parms{i}.Ref_Subsets;
       test_set{1} = cfd(i).stimfile;
    end
end

% Train on everything else by default
for i = 1:len;
    if ~isequal(cfd(i).stimfile, test_set{1})
        training_set{end+1} = cfd(i).stimfile;
    end
end

% TODO: Consider saving the parm/perf info to a global for later use?
% TODO: Consider moving parm/perf outside this fn?

log_dbg('Training sets selected: %s', ...
    strtrim(sprintf('%s ', training_set{:})));
log_dbg('Test set selected: %s', char(test_set{:}));

% ------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA SELECTION GUI

function query_db_button_Callback(hObject, eventdata, handles)
log_dbg('query_db_button pressed.');

% Clear the global variable storing state for the GUI (GS = 'Gui State')
global GS; GS = {};

% Invalidate data_selection_table, and other parts of the GUI
set(handles.data_selection_table, 'Data', {}); drawnow;

% Query the DB and select the train/test sets
GS.cellid = get(handles.cellid_text, 'String');
GS.cfd = query_db(GS.cellid);
[GS.training_set, GS.test_set] = select_train_test_sets(GS.cfd);

% Convert the above into a GUI viewable cell array
c = cell(length(GS.cfd), 3);
for i = 1:length(GS.cfd);
    c{i,1} = GS.cfd(i).stimfile; 
    c{i,2} = sum(ismember(GS.training_set, c{i,1})) > 0; 
    c{i,3} = sum(ismember(GS.test_set, c{i,1})) > 0;
end
set(handles.data_selection_table, 'Data', c); drawnow;

% TODO: Request that the files be loaded
% TODO: Load the first of the training sets
% TODO: Display the stimuli data.


% ------------------------------------------------------------------------
function update_raw_stim_plot(hObject, eventdata, handles)
global TIME STIM SAMPFREQ;
axes(handles.stim_view_axes); cla;
nsel = cellstr(get(handles.raw_stim_view_popup, 'String'));
plottype = nsel{get(handles.raw_stim_view_popup, 'Value')};
nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});
switch plottype
    case 'Time Series View'
        plot(TIME, STIM(stim_idx,:), 'k-');
        axis tight;
    case 'Spectrogram View'
        % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
        nwin = 4048;
        logfsgram(STIM(stim_idx,:)', nwin, SAMPFREQ, [], [], 500, 12);
        %caxis([-20,40]); % TODO: use a 'smarter' caxis here which discards
        % information based on a histogram to get rid of outliers.
end


function update_raw_resp_plot(hObject, eventdata, handles)
global TIME RESPAVG;
nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});
axes(handles.resp_view_axes); cla;
bar(TIME, RESPAVG(stim_idx,:), 0.01,'k');
axis tight;


function raster_freq_CreateFcn(hObject, eventdata, handles)
function raster_freq_Callback(hObject, eventdata, handles)
% TODO: Save the raster frequency to a global somewhere
%frq = str2double(get(hObject,'String'));

% 
% function select_training_set_button_Callback(hObject, eventdata, handles)
% global STIM RESP TIME RESPAVG SAMPFREQ;
% 
% index=1;
% options.includeprestim = 1;
% options.unit     = cfd(index).unit;
% options.channel  = cfd(index).channum;
% options.rasterfs = 100000; 
% SAMPFREQ = options.rasterfs; 
% 
% % TODO: make this load multiple stimulus and response files for a single
% % cellid
% 
% fprintf('Loading stimulus file: %s%s\n', cfd(index).stimpath, cfd(index).stimfile);
% stimfile = [cfd(index).stimpath cfd(index).stimfile];
% stim     = loadstimfrombaphy(stimfile, [], [], 'wav', options.rasterfs, 1, 0, options.includeprestim);
% 
% fprintf('Loading response file: %s/%s\n', cfd(index).path, cfd(index).respfile);
% respfile = [cfd(index).path cfd(index).respfile];
% [resp, tags] = loadspikeraster(respfile, options);
% 
% [d1 d2 d3] = size(stim);
% if d1 ~= 1
%     disp([d1 d2 d3]);
%     error('Stimulus matrix must initially have size 1xLxN, where L=length in samples, N=num of stimuli');
% end
% 
% STIM = permute(stim, [3 2 1]);  % Remove the irrelevant first dimension
% RESP = permute(resp, [3 1 2]);  % Rearrange to match (TODO: Try squeeze() instead)
% RESP(isnan(RESP)) = 0;          % Replace all NaN's with 0's. TODO: is this the right behavior?
% rsp = sum(RESP,3); 
% [e1 e2 e3] = size(rsp);
% RESPAVG = reshape(rsp, [e1 e2]); % Averaged response across all trials TODO: Squeeze()
% TIME = (1/options.rasterfs).*[1:d2]';
% 
% % TODO: Check that STIM, RESP, RESPAVG, TIME are all the proper sizes
% % PROBABLY this should be a function that you can call anytime to get
% % global values checked.
% 
% % Update other parts of the GUI
% set(handles.training_sets_listbox, 'String', char(cfd.stimfile));   
% set(handles.cellid_text, 'String', char(cellids(1)));   
% 
% % Update the popup for selecting a particular stimuli to graph
% c = {};
% for i = 1:e1;
%     c{i} = sprintf('%d',i);
% end
% set(handles.selected_stimuli_popup, 'String', char(c)); 
% set(handles.selected_stimuli_popup, 'Value', 1);
% 
% update_raw_stim_plot(hObject, eventdata, handles);
% update_raw_resp_plot(hObject, eventdata, handles);
%------------------------------------------------------------------------

function view_strf_button_Callback(hObject, eventdata, handles)
% TODO: If a TOR file exists in data_selection_table, use it to view
%strf_offline2();
% TODO: If a TOR file does NOT exist, print a dialog box saying so

function raw_stim_view_popup_CreateFcn(hObject, eventdata, handles)
function raw_stim_view_popup_Callback(hObject, eventdata, handles)
update_raw_stim_plot(hObject, eventdata, handles);

function raw_resp_view_popup_CreateFcn(hObject, eventdata, handles)
function raw_resp_view_popup_Callback(hObject, eventdata, handles)
update_raw_resp_plot(hObject, eventdata, handles);

function training_sets_listbox_CreateFcn(hObject, eventdata, handles)
function training_sets_listbox_Callback(hObject, eventdata, handles)
% TODO: change the graphs of everything!

function selected_stimuli_popup_CreateFcn(hObject, eventdata, handles)
function selected_stimuli_popup_Callback(hObject, eventdata, handles)
update_raw_stim_plot(hObject, eventdata, handles);
update_raw_resp_plot(hObject, eventdata, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPROCESSING WIDGETS

function update_prefilter_plots(handles)
global TIME STIM PF_COEFS RESPAVG PF_STIM SAMPFREQ;

nvals = cellstr(get(handles.preprocessing_view_popup, 'String'));
plottype = nvals{get(handles.preprocessing_view_popup, 'Value')};

nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});

n_filts = length(PF_COEFS);

axes(handles.preproc_view_axes); cla;
hold on;
switch plottype
    case 'Freq. Response'
        for filt_idx = 1:n_filts
            % Old way
            % [sos,g] = tf2sos(PF_COEFS{i}{1}, PF_COEFS{i}{2}); 
            % Hd = dfilt.df2tsos(sos,g);  % Create a dfilt object
            % h = fvtool(Hd);              % Plot magnitude response
            % set(h,'Analysis','magnitude') % Freq Magnitude response
            
            % New way: 
            ww = 0:(pi/1000):pi;
            H = freqz(PF_COEFS{filt_idx}{1}, PF_COEFS{filt_idx}{2}, ww);
            loglog(ww, abs(H), pickcolor(filt_idx));
            setAxisLabelCallback(gca, @(f) (f*SAMPFREQ/(3.14*2)), 'X');
            axis tight;
            % Best way?
            % P = bodeoptions;
            % P.PhaseVisible = 'off';
            % P.FreqUnits = 'Hz'; 
            % h = bodeplot(tf(PF_COEFS{filt_idx}{2}, PF_COEFS{filt_idx}{1}),P);               
        end
    case 'Filtered Stimulus'
        for filt_idx = 1:n_filts
            plot(TIME, squeeze(PF_STIM(stim_idx,:,filt_idx)), pickcolor(filt_idx));
            axis tight;
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

function refresh_preproc_view_button_Callback(hObject, eventdata, handles)
update_prefilter_plots(handles);

% TODO: Make this work for a ANY prefilter selection and move this to a
% prefilter model directory
function apply_prefilter_button_Callback(hObject, eventdata, handles)
global STIM SAMPFREQ PF_COEFS PF_STIM;
[d1, d2] = size(STIM);

% Get the low-pass and high-pass band limits
lp = eval(get(handles.low_freqs, 'String'));
hp = eval(get(handles.high_freqs, 'String'));
sf = eval(get(handles.smoothing_frqs, 'String'));

% TODO: Check that the inputs are valid
% if ~(isvector(lp) & isvector(hp) & isvector(sf))
%     set(handles.prefilter_status, 'String', 'ERROR: ALL PREFILTER SETTINGS MUST BE VECTORS AND OF THE SAME LENGTH!');
%     return;
% end

n_filts = length(lp);

% Make the bank of filters. I use a cell array here to allow different
% filter sizes and orders, if such a case actually comes up. 
PF_COEFS={}; 
for i = 1:n_filts
    [B,A] = ellip(4,0.5,50,[lp(i)/SAMPFREQ*2, hp(i)/SAMPFREQ*2]);   
    PF_COEFS{i} = {B,A};
end
% Filter the data
%set(handles.prefilter_status, 'String', 'Filtering...');
PF_STIM=[];
for i = 1:n_filts
    pfstmp1 = filter(PF_COEFS{i}{1}, PF_COEFS{i}{2}, STIM,[],2);
    pfstmp2 = abs(pfstmp1);
    % OPTIONAL LOW_PASS SMOOTHING ALGORITHM
    if sf ~= 0
        [B,A] = ellip(6,3,50,[sf(i)/SAMPFREQ*2]);
        pfs = filter(B,A, pfstmp2,[],2);
    else
        pfs = pfstmp2;
    end
    PF_STIM = cat(3, PF_STIM, pfs);
end

% Plot either the data or the frequency response
%set(handles.prefilter_status, 'String', 'Plotting...');
update_prefilter_plots(handles);
%set(handles.prefilter_status, 'String', 'Done.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CLASS WIDGETS
function update_fir_model(handles)
global PF_COEFS FIRHIST FIRBINSIZE FIRCOEFS;
FIRHIST = str2num(get(handles.fir_history,'String'));
FIRBINSIZE = str2num(get(handles.bin_size,'String'));
d = length(PF_COEFS); % Since this is a cell array, this returns just nfilts
l = ceil(FIRHIST/FIRBINSIZE);
FIRCOEFS = zeros(d,l);

function downsample_stimresp()
global PF_STIM RESPAVG SAMPFREQ DS_FREQ DS_STIM DS_TIME DS_RESPAVG FIRBINSIZE;
DS_FREQ = 1000 / FIRBINSIZE;  
disp('Downsampling the PF_STIM and RESP signals...');
DS_RESPAVG = conv_fn(RESPAVG, 2, @sum, SAMPFREQ/DS_FREQ, 0);
[d,l] = size(DS_RESPAVG);
DS_TIME = [0:1/DS_FREQ:l/DS_FREQ-1/DS_FREQ];
DS_STIM = conv_fn(PF_STIM, 2, @mean, SAMPFREQ/DS_FREQ, 0);



function make_predictions()
% Apply the FIR filters to corresponding downsampled stimuli to get the model prediction
% Since it is linear, the prediction is just the sum of the filters
% We assume that there are no second order terms combining elements of both filters

global DS_STIM DS_PREDS FIRCOEFS DS_PRED;
TOTAL_INPUT = sum(DS_STIM, 3);
DS_PREDS = [];

%disp('Applying FIR filters');
[n_filts, l] = size(FIRCOEFS);

for filt_idx = 1:n_filts 
    DS_PREDS = cat(3, DS_PREDS, sqrt(abs(filter(FIRCOEFS(filt_idx,:), 1, squeeze(DS_STIM(:,:,filt_idx)), [],2))));
end

DS_PRED = squeeze(sum(DS_PREDS, 3)); 





function plot_model(handles)
global FIRCOEFS DS_TIME DS_STIM DS_RESPAVG DS_PREDS DS_PRED FIRBINSIZE;

nvals = cellstr(get(handles.model_view_popup, 'String'));
plottype = nvals{get(handles.model_view_popup, 'Value')};

nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});

[n_filts, n_coefs] = size(FIRCOEFS);

axes(handles.model_view_axes); cla;
hold on;
switch plottype
    case 'FIR Shape'
        for filt_idx = 1:n_filts
            stem([1:n_coefs], FIRCOEFS(filt_idx,:), pickcolor(filt_idx));
        end
        setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;
    case 'Downsampled Stims'
        for filt_idx = 1:n_filts
            plot(DS_TIME, DS_STIM(stim_idx,:,filt_idx), pickcolor(filt_idx));
        end
        setAxisLabelCallback(gca, @(t) (t), 'X');
        axis tight; 
    case 'Downsampled Preds'
        make_predictions();
        for filt_idx = 1:n_filts
            plot(DS_TIME, DS_PREDS(stim_idx,:,filt_idx), pickcolor(filt_idx));
        end
        setAxisLabelCallback(gca, @(t) (t), 'X');
        axis tight;
    case 'Downsampled Resp'
        %plot(DS_TIME, DS_RESPAVG(stim_idx,:));
        s1 = 1/max(DS_RESPAVG(stim_idx,:));
        s2 = 1/max(DS_PRED(stim_idx,:));
            plot(DS_TIME, s1*DS_RESPAVG(stim_idx,:), pickcolor(0), ...
             DS_TIME, s2*DS_PRED(stim_idx,:), pickcolor(1));
%        plot(DS_TIME, s1*DS_RESPAVG(stim_idx,:), pickcolor(0));
        setAxisLabelCallback(gca, @(t) (t), 'X');
        axis tight; 
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

function downsample_button_Callback(hObject, eventdata, handles)
update_fir_model(handles);
downsample_stimresp();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMIZATION WIDGETS

function x = pak_FIRCOEFS(FC)
x = FC(:)';

function FC = unpak_FIRCOEFS(x)
global PF_COEFS FIRHIST FIRBINSIZE;
d = length(PF_COEFS);
l = ceil(FIRHIST/FIRBINSIZE);
FC = reshape(x, d, l);


function z = correlation_of_downsampled_signals(x)
% Returns the correlation of the binned stimulus and binned response
global DS_RESPAVG DS_PREDS DS_PRED FIRCOEFS;
FIRCOEFS = unpak_FIRCOEFS(x);  % Update the global
make_predictions(); % TODO: This should be more functional!

[n_stims, n_ds_samps] = size(DS_RESPAVG);

corrs = zeros(1, n_stims);

% The average correlation across ALL trials is what we care about
% for stim_idx = 1:n_stims
%     R = corrcoef(DS_PRED(stim_idx,:), DS_RESPAVG(stim_idx,:));
%     R(isnan(R)) = 0; % Replace NaN's with 0 correlations again.
%     corrs(stim_idx) = R(2,1);
% end
%z = mean(corrs);

% TODO: Allow option to do correlation average or concatenated.
V1 = reshape(DS_PRED.',[],1);
V2 = reshape(DS_RESPAVG.',[],1);
R = corrcoef(V1,V2);
R(isnan(R)) = 0;
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
global FIRCOEFS;
% Get the number of iterations
n_iters = str2num(get(handles.termination_iterations, 'String'));

% Get the starting score of the current filter
x_0 = pak_FIRCOEFS(FIRCOEFS);

stepsize = 1.0;
[x_bst, s_bst] = boosting(x_0, @correlation_of_downsampled_signals, @(n,x,s)(n > n_iters), stepsize);
unpak_FIRCOEFS(x_bst);


% ------------------------
function initialize_optplot_menu(handle)
menuopts = {'Param Likelihoods', 'Pred/PSTH scatter', ...
        'KS Plot', 'Raw ISI', 'Time-Scaled ISI', ...
        'Cond. Intensity Fn', 'Time Scaling Fn'};

set(handle, 'String', char(menuopts));

function setoptplot(handles, plothandle, plottype)
global DS_RESPAVG DS_PRED;

nvals = cellstr(get(handles.selected_stimuli_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stimuli_popup, 'Value')});

axes(plothandle); cla;
hold on;
switch plottype
    case 'Param Likelihoods'
        %plot(DS_PRED, DS_RESPAVG, 'k.');
        %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;
    case 'Pred/PSTH scatter'
        plot(DS_PRED(stim_idx,:), DS_RESPAVG(stim_idx,:), 'k.');
        %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;
    case 'KS Plot'
        %plot(DS_PRED, DS_RESPAVG, 'k.');
        %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;
    case 'Raw ISI'
        %plot(DS_PRED, DS_RESPAVG, 'k.');
        %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;        
    case 'Time-Scaled ISI'
        %plot(DS_PRED, DS_RESPAVG, 'k.');
        %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;
    case 'Cond. Intensity Fn'
        %plot(DS_PRED, DS_RESPAVG, 'k.');
        %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;        
    case 'Time Scaling Fn'
        %plot(DS_PRED, DS_RESPAVG, 'k.');
        %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
        axis tight;
end
hold off;


function optplot1popup_CreateFcn(hObject, eventdata, handles)
initialize_optplot_menu(hObject);
function optplot1popup_Callback(hObject, eventdata, handles)
nvals = cellstr(get(handles.optplot1popup, 'String'));
plottype = nvals{get(handles.optplot1popup, 'Value')};
setoptplot(handles, handles.optplot1, plottype);

function optplot2popup_CreateFcn(hObject, eventdata, handles)
initialize_optplot_menu(hObject);
function optplot2popup_Callback(hObject, eventdata, handles)
nvals = cellstr(get(handles.optplot2popup, 'String'));
plottype = nvals{get(handles.optplot2popup, 'Value')};
setoptplot(handles, handles.optplot2, plottype);

function optplot3popup_CreateFcn(hObject, eventdata, handles)
initialize_optplot_menu(hObject);
function optplot3popup_Callback(hObject, eventdata, handles)
nvals = cellstr(get(handles.optplot3popup, 'String'));
plottype = nvals{get(handles.optplot3popup, 'Value')};
setoptplot(handles, handles.optplot3, plottype);

function optplot4popup_CreateFcn(hObject, eventdata, handles)
initialize_optplot_menu(hObject);
function optplot4popup_Callback(hObject, eventdata, handles)
nvals = cellstr(get(handles.optplot4popup, 'String'));
plottype = nvals{get(handles.optplot4popup, 'Value')};
setoptplot(handles, handles.optplot4, plottype);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REPORTING WIDGETS

% Nothing so far!
function select_test_set_button_Callback(hObject, eventdata, handles)

function dump_data_button_Callback(hObject, eventdata, handles)

tempdata = cell(4,3);
tempdata{1,1} = true;
tempdata{2,1} = false;
tempdata{1,2} = 'Bandpass Lo Frq'; tempdata{1,3} = 9000;
tempdata{2,2} = 'Bandpass Hi Frq'; tempdata{2,3} = 14000;

set(handles.preproc_data_table, 'Data', tempdata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on selection change in popupmenu17.
function popupmenu17_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu17 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu17


% --- Executes during object creation, after setting all properties.
function popupmenu17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
