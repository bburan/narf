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

% Last Modified by GUIDE v2.5 21-Nov-2012 14:08:42

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

% Create a scrollable edit box in the status panel
hEdit = uicontrol(handles.uipanel6, 'Style','edit', 'FontSize',9, ...
    'Min',0, 'Max',2, 'HorizontalAlignment','left', ...
    'Units','normalized', 'Position',[0 0 1 1], ...
    'String','GUI Initialized');

%%%%%%%%%%%%%% GLOBAL STUFF BELOW %%%%%%%%%%%%%%%%%%

% LOG FILE RELATED STUFF
global LOG_HANDLE LOG_LENGTH LOG_BUFFER LOG_LEVEL;
% TODO: Make these more editable at the top level?
LOG_LEVEL = 0;       % 0=debug, 1=informative, 2=normal, 3=warnings, 4=err 
LOG_HANDLE = hEdit;  % Where to print log strings
LOG_LENGTH = 6;      % 6 lines visible
LOG_BUFFER = cell(1,LOG_LENGTH);
for i = 1:LOG_LENGTH
    LOG_BUFFER{i} = '';
end

% GLOBAL GUI STATE (GS) STRUCTURE 
global GS; 
GS = [];

% PATHS
global NARF_PATH PREPROC_DIR DOWNSAMP_DIR MODEL_DIR ...
    STOCHAST_DIR SAMPLING_DIR PERF_METRIC_DIR TERMINATION_DIR;
NARF_PATH       = '/home/ivar/matlab/narf/';
MODULES_DIR     = 'modules/'
SAMPLING_DIR    = 'optim_0_sampling/';
PERF_METRIC_DIR = 'optim_1_perf_metric/';
TERMINATION_DIR = 'optim_2_termination/';

%%%%%%%%%%%%%% GLOBAL STUFF ABOVE %%%%%%%%%%%%%%%%%%

% Add necessary directories to NARF's path
addpath([NARF_PATH filesep 'utils'], ...
        [NARF_PATH filesep MODULES_DIR], ...
        [NARF_PATH filesep SAMPLING_DIR], ...
        [NARF_PATH filesep PERF_METRIC_DIR], ...
        [NARF_PATH filesep TERMINATION_DIR]);

% Build a list of the files in the 'modules/' directory
mods = scan_directory_for_modules('~/matlab/narf/modules/');

% Add the model pane
narf_modelpane(handles.model_structure_panel, mods);
    
% % Invalidate all data tables
set(handles.data_selection_table, 'Data', {});
set(handles.sampling_data_table, 'Data', {});
set(handles.perf_metric_data_table, 'Data', {});
set(handles.term_cond_data_table, 'Data', {});
drawnow;

% --- Outputs from this function are returned to the command line.
function varargout = narf_gui_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USEFUL LOGGING FUNCTIONS

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
function cfd = query_db(cellid)
log_dbg('query_db(''%s'');', cellid);

% Returns a list of raw files with this cell ID
[cfd, cellids, cellfileids] = dbgetscellfile('cellid', cellid);

% If there is not exactly one cell file returned, throw an error.
if ~isequal(length(cellids), 1)
    log_err('BAPHY gave %d cellids yet I want only 1.', length(cellids));
end

% ------------------------------------------------------------------------
function [training_set, test_set] = select_train_test_sets(cfd)
log_dbg('Selecting training and test sets');
len = length(cfd);
training_set = {};
test_set = {};
test_set_reps = 0;
parms = cell(1, len);
perfs = cell(1, len);

% Load parms, perfs. Select set with the most repetitions as test set.
% TODO: Right now this is not really choosing the most repetitions, just
% the number of Ref_Subsets, which just works by coincidence?
for i = 1:len;
    [parms{i}, perfs{i}] = dbReadData(cfd(i).rawid);
    if (isfield(parms{i}, 'Ref_Subsets') & ...
            parms{i}.Ref_Subsets > test_set_reps)
       test_set_reps = parms{i}.Ref_Subsets;
       test_set{1} = cfd(i).stimfile;
    end
end

% Train on every other passive trial by default
for i = 1:len;
    if (~isequal(cfd(i).stimfile, test_set{1}) & ...
         isequal(cfd(i).behavior, 'passive'))
        training_set{end+1} = cfd(i).stimfile;
    end
end

% TODO: Consider saving the parm/perf info to a global for later use?
% TODO: Consider moving parm/perf outside this fn?
log_dbg('Training sets selected: %s', ...
    strtrim(sprintf('%s ', training_set{:})));
log_dbg('Test set selected: %s', char(test_set{:}));

%------------------------------------------------------------------------
function generic_model_data_table_update(hObject, GSfield, GSselected_name)
% Since the data tables look the same for preprocessing, downsampling,
% model, and stochasticity, let's abstract their similarities here.
% Does two things:
% 1. Pull out the present values and update the params struct first
% 2. Update which parameters are desired to be fit with the optimization
global GS;
s = extract_field_val_pairs(hObject, 2, 3);
fns = fieldnames(s);
dt = GS.(GSfield).(GS.(GSselected_name));
for i = 1:length(fns);
    dt.params.(fns{i}) = s.(fns{i});
end
GS.(GSfield).(GS.(GSselected_name)).params = dt.fn(dt.params); 

% 2. Update which parameters are desired to be fit with the optimization
GS.(GSfield).(GS.(GSselected_name)).fittable_params = ...
     extract_checked_fields(hObject, 1, 2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA SELECTION GUI

function query_db_button_Callback(hObject, eventdata, handles)
log_dbg('query_db_button pressed.');

% GS: Clear various parts of the global struct
global GS;

% Invalidate data_selection_table, and other parts of the GUI
set(handles.data_selection_table, 'Data', {}); drawnow;
% axes(handles.stim_view_axes); cla;
% axes(handles.resp_view_axes); cla;
% axes(handles.preproc_view_axes); cla;
% axes(handles.downsamp_view_axes); cla;
% axes(handles.optplot1); cla;
% axes(handles.optplot2); cla;
% axes(handles.optplot3); cla;
% set(handles.selected_stimfile_popup, 'String', '');   
% set(handles.selected_stim_idx_popup, 'String', '');

% Query the DB and select the train/test sets
GS.cellid = get(handles.cellid_text, 'String');
GS.cfd = query_db(GS.cellid);
[GS.training_set, GS.test_set] = select_train_test_sets(GS.cfd);

% Convert the above into a GUI viewable cell array, then update GUI
c = cell(length(GS.cfd), 3);
for i = 1:length(GS.cfd);
    c{i,1} = GS.cfd(i).stimfile; 
    c{i,2} = sum(ismember(GS.training_set, c{i,1})) > 0; 
    c{i,3} = sum(ismember(GS.test_set, c{i,1})) > 0;
end
set(handles.data_selection_table, 'Data', c); drawnow;

% ------------------------------------------------------------------------
function data_selection_table_CellEditCallback(hObject, eventdata, handles)
global GS;
log_dbg('data_selection_table was modified, updating GB.*_set.');

% Refresh GB.training_set and GB.test_set from the GUI 
d = get(hObject, 'Data');
[r, c] = size(d);
training_set = {};
test_set = {};
j = 1;
k = 1;
for i = 1:r
    if d{i,2},
        training_set{j} = d{i,1};
        j = j+1;
    elseif d{i,3}
        test_set{k} = d{i,1};
        k = k+1;
    end
end

GS.training_set = training_set;
GS.test_set = test_set;

% ------------------------------------------------------------------------
function load_above_files_button_Callback(hObject, eventdata, handles)
global GS;
set(handles.selected_stimfile_popup, 'String', '');  
set(handles.selected_stim_idx_popup, 'String', '');
drawnow;
GS.dat = load_stim_resps(GS.cfd, GS.training_set, GS.test_set);
update_selected_stimfile_popup(handles); % Push changes to GUI
update_selected_stim_idx_popup(handles);
selected_stimfile_popup_Callback([], [], handles); % Trigger GUI callbacks
selected_stim_idx_popup_Callback([], [], handles);
raw_stim_view_popup_Callback([],[],handles); % Will update
raw_resp_view_popup_Callback([],[],handles);

% ------------------------------------------------------------------------
function update_selected_stimfile_popup(handles)
global GS;
fns = fieldnames(GS.dat);
set(handles.selected_stimfile_popup, 'String', char(fns));   
GS.selected_stimfile = fns{1};

% ------------------------------------------------------------------------
function update_selected_stim_idx_popup(handles)
global GS;
stimfile = GS.selected_stimfile;
c = {}; 
if isfield(GS.dat, stimfile)
    [d1, d2] = size(GS.dat.(stimfile).raw_stim);
    for i = 1:d1
        c{i} = sprintf('%d',i); 
    end
    set(handles.selected_stim_idx_popup, 'String', char(c)); 
    set(handles.selected_stim_idx_popup, 'Value', 1);
else
    log_err('Selected stimulus file not found: %s', stimfile);
end

% ------------------------------------------------------------------------
function update_raw_resp_plot(handles)
log_dbg('Updating raw_resp_plot');
global GS;
 
% Only update if all fields are defined
if ~(isfield(GS, 'raw_resp_plot_type') & ...
     isfield(GS, 'selected_stim_idx') & ...
     isfield(GS, 'selected_stimfile'))
    log_dbg('Ignoring raw resp plot update since not all fields ready.');
    return;
end

plottype = GS.raw_resp_plot_type;
stim_idx = GS.selected_stim_idx;
stimfile = GS.selected_stimfile;
obj = GS.dat.(stimfile);
[S, N, R] = size(obj.raw_resp);

axes(handles.resp_view_axes); cla; hold on;
switch plottype
    case 'Time Series View'
        [xs,ys] = find(obj.raw_respavg(stim_idx, :) > 0);
        bar(obj.raw_time(ys), obj.raw_respavg(stim_idx,ys), 0.01, 'k');
        axis([0 obj.raw_time(end) 0 2])
    case 'Dot Rasterization'
        for j = 1:R
            [xs,ys] = find(obj.raw_resp(stim_idx, :, j) > 0);
            plot(obj.raw_time(ys), j/R*obj.raw_resp(stim_idx,ys,j), 'k.');
        end
        axis([0 obj.raw_time(end) 1/(2*R) 1+(1/(2*R))])
end
hold off;

%------------------------------------------------------------------------
function selected_stimfile_popup_CreateFcn(hObject, eventdata, handles)
function selected_stimfile_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(handles.selected_stimfile_popup, 'String'));
stimfile = nvals{get(handles.selected_stimfile_popup, 'Value')};
GS.selected_stimfile = stimfile;
update_selected_stim_idx_popup(handles);
update_raw_stim_plot(handles);
update_raw_resp_plot(handles);
preproc_view_popup_Callback(handles.preproc_view_popup, [], handles);
update_model_view_plot(handles);

%------------------------------------------------------------------------
function selected_stim_idx_popup_CreateFcn(hObject, eventdata, handles)
function selected_stim_idx_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(handles.selected_stim_idx_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stim_idx_popup, 'Value')});
GS.selected_stim_idx = stim_idx;
update_raw_stim_plot(handles);
update_raw_resp_plot(handles);
preproc_view_popup_Callback(handles.preproc_view_popup, [], handles);
update_model_view_plot(handles);

%------------------------------------------------------------------------
function raw_stim_view_popup_CreateFcn(hObject, eventdata, handles)
function raw_stim_view_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(handles.raw_stim_view_popup, 'String'));
plottype = nvals{get(handles.raw_stim_view_popup, 'Value')};
GS.raw_stim_plot_type = plottype;
update_raw_stim_plot(handles);

%------------------------------------------------------------------------
function raw_resp_view_popup_CreateFcn(hObject, eventdata, handles)
function raw_resp_view_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(handles.raw_resp_view_popup, 'String'));
plottype = nvals{get(handles.raw_resp_view_popup, 'Value')};
GS.raw_resp_plot_type = plottype;
update_raw_resp_plot(handles);

%------------------------------------------------------------------------
function view_strf_button_Callback(hObject, eventdata, handles)
log_dbg('view_strf_button pressed');
global GS;
for i = 1:length(GS.cfd);
    % TODO: Replace magic number  1 with better description of TOR files
    if (GS.cfd(i).runclassid == 1) 
        log_dbg('Making STRF for ''%s''', GS.cfd(i).stimfile);
        figure; 
        strf_offline2([GS.cfd(i).stimpath GS.cfd(i).stimfile], ...
            [GS.cfd(i).path GS.cfd(i).respfile], ...
            GS.cfd(i).channum, GS.cfd(i).unit);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPROCESSING WIDGETS

%------------------------------------------------------------------------
function preproc_popup_CreateFcn(hObject, eventdata, handles)
function preproc_popup_Callback(hObject, eventdata, handles)
generic_model_selecting_popup(hObject, ...
    'preproc', 'selected_preproc_name', handles.preproc_data_table);

%------------------------------------------------------------------------
function preproc_data_table_CellEditCallback(hObject, eventdata, handles)
global GS;
log_dbg('preproc_data_table modified. Click ''preprocess!'' to refresh.');
generic_model_data_table_update(hObject, 'preproc', 'selected_preproc_name');
GS.dat.(GS.selected_stimfile).pp_stim = []; % Invalidate data
axes(handles.preproc_view_axes); cla;       % Clear plot

%------------------------------------------------------------------------
function preproc_button_Callback(hObject, eventdata, handles)
global GS;
% TODO: Check that data is already loaded and is available for preprocing!
% TODO: Check that the object dhas a .preproc_fn method!
log_inf('Preprocessing...');
%-------
% TODO: Move this core computation somewhere else!
f = fieldnames(GS.dat); 
spp = GS.preproc.(GS.selected_preproc_name).params;
for i = 1:length(f)
    % Feed data into that new filter
    GS.dat.(f{i}).pp_stim = spp.preproc_fn(spp, GS.dat.(f{i}).raw_stim); 
end
% TODO: Check that the length of the preprocessed vector is the SAME size
% as the original raw stimulus.
%-------
log_inf('Done preprocessing.');
% Reset the GUI using callbacks
set(handles.preproc_index_popup, 'String', 'NONE');
preproc_index_popup_Callback(handles.preproc_index_popup, [], handles);
preproc_view_popup_Callback(handles.preproc_view_popup, [], handles);

%------------------------------------------------------------------------
function preproc_view_popup_CreateFcn(hObject, eventdata, handles)
function preproc_view_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(hObject, 'String'));
plottype = nvals{get(hObject, 'Value')};
GS.preproc_view_plot_type = plottype;
update_preproc_view_plot(handles);
update_downsamp_view_plot(handles);

%------------------------------------------------------------------------                                
function preproc_index_popup_CreateFcn(hObject, eventdata, handles)
function preproc_index_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(hObject, 'String'));
GS.selected_preproc_idx = str2num(nvals{get(hObject, 'Value')});
if isempty(GS.selected_preproc_idx)
    GS.selected_preproc_idx = 1;
end
update_preproc_view_plot(handles);
update_downsamp_view_plot(handles);

%------------------------------------------------------------------------
function update_preproc_view_plot(handles)
log_dbg('Updating preproc_view_plot');
global GS;

% Only update if all fields are defined
if ~(isfield(GS, 'selected_stimfile') & ...
     isfield(GS, 'selected_preproc_name') & ...
     isfield(GS, 'selected_preproc_idx') & ...    
     isfield(GS, 'preproc_view_plot_type') & ...
     isfield(GS, 'dat') & ...
     isfield(GS.dat, GS.selected_stimfile) & ...
     isfield(GS.dat.(GS.selected_stimfile), 'pp_stim'))
   log_dbg('Ignoring preproc plot update since not all fields ready.');
   return  
end

dat = GS.dat.(GS.selected_stimfile);
[n_stims, n_samps, n_filts] = size(dat.pp_stim);

% By default, set the preproc_index_popup to be disabled but populated
c={};
for i = 1:n_filts
    c{i} = sprintf('%d', i);
end    
set(handles.preproc_index_popup, 'String', char(c));
set(handles.preproc_index_popup, 'Enable', 'Off');
set(handles.view_preproc_index_label, 'Enable', 'Off');

spp = GS.preproc.(GS.selected_preproc_name);

axes(handles.preproc_view_axes); cla;
hold on;
switch GS.preproc_view_plot_type
    case 'Frequency Response'
        % If the filter has a frequency response method defined, call it
        if isfield(spp.params, 'freq_resp_plot_fn')
            spp.params.freq_resp_plot_fn(spp.params);
        else
            log_dbg('No fn found to plot freq response.');
        end
    case 'Filtered Stimulus'
        for filt_idx = 1:n_filts
             plot(dat.raw_time, ...
                  squeeze(dat.pp_stim(GS.selected_stim_idx,:,filt_idx)), ...
                  pickcolor(filt_idx));
        end
        setAxisLabelCallback(gca, @(t) (t), 'X');
        axis tight;
     case 'Filtered Spectrogram'
        set(handles.preproc_index_popup, 'Enable', 'On');
        set(handles.view_preproc_index_label, 'Enable', 'On');
        logfsgram(dat.pp_stim(GS.selected_stim_idx,:,GS.selected_preproc_idx)', 4048, 100000, [], [], 500, 12); 
        caxis([-20,40]);
        setAxisLabelCallback(gca, @(t) (t), 'X');
        drawnow;
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOWNSAMPLING WIDGETS

%------------------------------------------------------------------------
function downsamp_popup_CreateFcn(hObject, eventdata, handles)
function downsamp_popup_Callback(hObject, eventdata, handles)
generic_model_selecting_popup(hObject, ...
    'downsamp', 'selected_downsamp_name', handles.downsamp_data_table);

%------------------------------------------------------------------------
function downsamp_data_table_CellEditCallback(hObject, eventdata, handles)
global GS;
log_dbg('downsampling_data_table modified. Click ''downsample!'' to refresh plots.');
generic_model_data_table_update(hObject, 'downsamp', 'selected_downsamp_name');
GS.dat.(GS.selected_stimfile).ds_stim = []; % Invalidate data
axes(handles.downsamp_view_axes); cla;    % Clear plot

%------------------------------------------------------------------------
function downsample_button_Callback(hObject, eventdata, handles)
global GS;
% TODO: Check that all data is ready for the downsampling!
% TODO: Check that the object has a .downsamp_fn method!
log_inf('Downsampling...');
% TODO: Move this core computation somewhere else!
f = fieldnames(GS.dat); 
p = GS.downsamp.(GS.selected_downsamp_name).params;
fs = p.raw_freq;
dsfs = p.ds_freq;
ratio = ceil(fs/dsfs);
% Process the stim, resp, and time matrices for each file
for i = 1:length(f)
    [ds_stim, ds_respavg] = p.downsamp_fn(p, GS.dat.(f{i}).pp_stim, ...
                                             GS.dat.(f{i}).raw_respavg); 
    GS.dat.(f{i}).ds_stim = ds_stim;
    GS.dat.(f{i}).ds_respavg = ds_respavg; 
    % TODO: How should the time resampling occur?
    GS.dat.(f{i}).ds_time = linspace(1/dsfs, GS.dat.(f{i}).raw_time(end), ...
                                     length(GS.dat.(f{i}).ds_respavg));
end
log_inf('Done downsampling.');
% Reset the GUI using callbacks
downsamp_view_popup_Callback(handles.downsamp_view_popup, [], handles);

%------------------------------------------------------------------------
function downsamp_view_popup_CreateFcn(hObject, eventdata, handles)
function downsamp_view_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(hObject, 'String'));
plottype = nvals{get(hObject, 'Value')};
GS.downsamp_view_plot_type = plottype;
update_downsamp_view_plot(handles);

%------------------------------------------------------------------------
function update_downsamp_view_plot(handles)
log_dbg('Updating downsamp_view_plot');
global GS;

% Only update if all fields are defined
if ~(isfield(GS, 'selected_stimfile') & ...
     isfield(GS, 'selected_downsamp_name') & ...
     isfield(GS, 'downsamp_view_plot_type') & ...
     isfield(GS, 'dat') & ...
     isfield(GS.dat, GS.selected_stimfile) & ...
     isfield(GS.dat.(GS.selected_stimfile), 'pp_stim') & ...
     isfield(GS.dat.(GS.selected_stimfile), 'ds_stim') & ...
     isfield(GS.dat.(GS.selected_stimfile), 'ds_respavg'))
    log_dbg('Ignoring downsample plot update since not all fields ready.');
   return  
end

% Get at the selected downsamp data
dat = GS.dat.(GS.selected_stimfile);
[n_stims, n_samps, n_filts] = size(dat.ds_stim);
sds = GS.downsamp.(GS.selected_downsamp_name);

% Depending on the type of downsamp plot desired
axes(handles.downsamp_view_axes); cla;
hold on;
switch GS.downsamp_view_plot_type
    case 'Downsampled Stimulus'
        for filt_idx = 1:n_filts
             plot(dat.ds_time, ...
                  squeeze(dat.ds_stim(GS.selected_stim_idx,:,filt_idx)), ...
                  pickcolor(filt_idx));
            axis tight;
        end
    case 'Downsampled Response'
        for filt_idx = 1:n_filts
             plot(dat.ds_time, ...
                  squeeze(dat.ds_respavg(GS.selected_stim_idx,:)), ...
                  pickcolor(filt_idx));
            axis tight;
        end
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CLASS WIDGETS

%------------------------------------------------------------------------
function model_popup_CreateFcn(hObject, eventdata, handles)
function model_popup_Callback(hObject, eventdata, handles)
generic_model_selecting_popup(hObject, ...
    'model', 'selected_model_name', handles.model_data_table);

%------------------------------------------------------------------------
function model_data_table_CellEditCallback(hObject, eventdata, handles)
global GS;
log_dbg('model_data_table modified. Click ''predict!'' to refresh plots.');
generic_model_data_table_update(hObject, 'model', 'selected_model_name');
GS.dat.(GS.selected_stimfile).ds_pred = []; % Invalidate data
axes(handles.model_view_axes); cla;    % Clear plot

%------------------------------------------------------------------------    
function model_button_Callback(hObject, eventdata, handles)
global GS;
% TODO: Check that all data is ready for the model fn to run!
% TODO: Check that the object has a .model_fn method!
log_inf('Running model...');
% TODO: Move this core computation somewhere else!
f = fieldnames(GS.dat); 
p = GS.model.(GS.selected_model_name).params;
% Process the stim, resp, and time matrices for each file
for i = 1:length(f)
    GS.dat.(f{i}).ds_pred = p.model_fn(p, GS.dat.(f{i}).ds_stim); 
end
log_inf('Done running model.');
% Reset the GUI using callbacks
model_view_popup_Callback(handles.model_view_popup, [], handles);


%------------------------------------------------------------------------
function model_view_popup_CreateFcn(hObject, eventdata, handles)
function model_view_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(hObject, 'String'));
plottype = nvals{get(hObject, 'Value')};
GS.model_view_plot_type = plottype;
update_model_view_plot(handles);

%------------------------------------------------------------------------
function update_model_view_plot(handles)
log_dbg('Updating model_view_plot');
global GS;

% Only update if all fields are defined
if ~(isfield(GS, 'selected_stimfile') & ...
     isfield(GS, 'selected_model_name') & ...
     isfield(GS, 'model_view_plot_type') & ...
     isfield(GS, 'dat') & ...
     isfield(GS.dat, GS.selected_stimfile) & ...
     isfield(GS.dat.(GS.selected_stimfile), 'ds_stim') & ...
     isfield(GS.dat.(GS.selected_stimfile), 'ds_pred') & ...
     isfield(GS.dat.(GS.selected_stimfile), 'ds_respavg'))
    log_dbg('Ignoring model plot update since not all fields ready.');
   return  
end

mod = GS.model.(GS.selected_model_name);
dat = GS.dat.(GS.selected_stimfile);
plottype = GS.model_view_plot_type;

axes(handles.model_view_axes); cla;
hold on;
switch plottype
	case 'Prediction vs. Reality'
        % Scale the response and prediction in case they have wildly
        % different scales (a common problem when using a correlation
        % coefficient-type performance metric is used to fit the model
        respavg = squeeze(dat.ds_respavg(GS.selected_stim_idx,:));
        rs = mean(respavg);
        stim = squeeze(dat.ds_pred(GS.selected_stim_idx,:));
        ss = mean(stim);
        % Plot 
        plot(dat.ds_time, (1/rs)*respavg, 'k-', ...
             dat.ds_time, (1/ss)*stim, 'r-');
        %setAxisLabelCallback(gca, @(t) (t), 'X');
        axis tight;
    case 'FIR Coefs'
        if isfield(mod.params, 'plot_coefs_fn')
            mod.params.plot_coefs_fn(mod.params);
        else
            log_dbg('No fn found to plot coefficients');
        end
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STOCHASTICITY AND ITS PLOTS


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

function stochast_popup_CreateFcn(hObject, eventdata, handles)
function stochast_popup_Callback(hObject, eventdata, handles)

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
menuopts = {'Pred/Resp Scatter', 'Raw ISI', 'Intensity-Scaled ISI', ...
        'KS Plot'};

set(handle, 'String', char(menuopts));

function setoptplot(handles, plothandle, plottype)
global GS;

stim_idx = GS.selected_stim_idx;
dat = GS.dat.(GS.selected_stimfile);

axes(plothandle); cla;
hold on;
switch plottype
    case 'Pred/Resp Scatter'
        plot(dat.ds_respavg(stim_idx,:), dat.ds_pred(stim_idx,:), 'k.');
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
    case 'KS Plot'
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
function check_stochast_button_Callback(hObject, eventdata, handles)

% --- Executes on button press in auto_recalc_checkbox.
function auto_recalc_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to auto_recalc_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_recalc_checkbox


% --- Executes on button press in save_model_params_button.
function save_model_params_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_model_params_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in load_model_params_button.
function load_model_params_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_model_params_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


