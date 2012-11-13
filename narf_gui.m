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

% Last Modified by GUIDE v2.5 12-Nov-2012 16:11:05

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
global NARF_PATH PREPROCESSING_DIR DOWNSAMPLING_DIR MODEL_DIR ...
    STOCHASTICITY_DIR SAMPLING_DIR PERF_METRIC_DIR TERMINATION_DIR;
NARF_PATH = '/home/ivar/matlab/narf/';
PREPROCESSING_DIR = 'stage_0_preprocessing/';
DOWNSAMPLING_DIR  = 'stage_1_downsampling/';
MODEL_DIR         = 'stage_2_model/';
STOCHASTICITY_DIR = 'stage_3_stochasticity/';
SAMPLING_DIR      = 'optim_0_sampling';
PERF_METRIC_DIR   = 'optim_1_perf_metric/';
TERMINATION_DIR   = 'optim_2_termination/';

%%%%%%%%%%%%%% GLOBAL STUFF ABOVE %%%%%%%%%%%%%%%%%%

% Add necessary directories to NARF's path
addpath([NARF_PATH filesep 'utils'], ...
        [NARF_PATH filesep PREPROCESSING_DIR], ...
        [NARF_PATH filesep DOWNSAMPLING_DIR], ...
        [NARF_PATH filesep MODEL_DIR], ...
        [NARF_PATH filesep STOCHASTICITY_DIR], ...
        [NARF_PATH filesep SAMPLING_DIR], ...
        [NARF_PATH filesep PERF_METRIC_DIR], ...
        [NARF_PATH filesep TERMINATION_DIR]);

% Invalidate all data tables
set(handles.data_selection_table, 'Data', {});
set(handles.sampling_data_table, 'Data', {});
set(handles.perf_metric_data_table, 'Data', {});
set(handles.term_cond_data_table, 'Data', {});
set(handles.preproc_data_table, 'Data', {});
set(handles.model_data_table, 'Data', {});
set(handles.stochasticity_data_table, 'Data', {});
set(handles.downsampling_data_table, 'Data', {});

% Initialize the preprocessing popup menu and data table
scanpath = fullfile(NARF_PATH, PREPROCESSING_DIR);
GS.preprocs = scan_directory_for_functions(scanpath);
fs = fieldnames(GS.preprocs);
GS.selected_preproc_name = fs{1}; % Use the first as the default
set(handles.preproc_popup, 'String', ...
    char(GS.preprocs.(fs{1}).params.pretty_name));
set(handles.preproc_popup, 'Value', 1);
set(handles.preproc_index_popup, 'String', '1'); % Initialize display idx too
set(handles.preproc_index_popup, 'Value', 1);
preproc_popup_Callback(handles.preproc_popup, [], handles);
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
% Core functionality of the system should be independent of the GUI

function run_narf(cellid, dataset_sel_fn, opt_fn, perf_metric, term_cond)
cfd = query_db(cellid);
[train_set, test_set] = select_train_test_sets(cfd); 
dat = load_stim_resps(cfd, train_set, test_set);
% Load the preprocessing model
% Load the preprocessing model params
%    OR Autodetect the bands of the SPN (when doing single-filter analysis)
% mod = load_model(preproc, model, stochasticity); % Load the model form
% Initialize the model OR
%    OR load the model parameters (save under cellid+stimfile in a 'saved' dir)

% preproc the data
% Downsample the filtered data
% Run the optimization on the downsampled data
% Return the best fitting parameters and/or save them

% ------------------------------------------------------------------------
function cfd = query_db(cellid)
log_dbg('query_db(''%s'');', cellid);

% Returns a list of raw files with this cell ID
[cfd, cellids, cellfileids] = dbgetscellfile('cellid', cellid);

% If there is not exactly one cell file returned, throw an error.
if ~isequal(length(cellids), 1)
    log_err('BAPHY gave me %d cellids yet I need one.', length(cellids));
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

% ------------------------------------------------------------------------
function dat = load_stim_resps(cfd, training_set, test_set)
% Return a 'dat' struct with 100Khz rasterized stims and responses. 'dat'
% struct has field names that are the stimfiles found in training_set and
% test_set, and so may be accessed using 'dat.mystimfilename.raw_stim', for
% example. The following fields will eventually be defined for the
% sub-structures held under each stimfile key:
%
%    raw_stim    The 100Khz rasterized stimulus signal    [SxN]
%    raw_resp    The 100Khz rasterized response signal    [SxNxR]
%    raw_respavg The 100Khz average response signal       [SxN]
%    raw_time    The 100Khz rasterized time index         [1xN]
%    pp_stim     The preprocessed stimulus                [SxNxF]
%    pp_resp     The preprocessed response                [SxNxRxF] 
%    pp_respavg  The preprocessed average response        [SxNxF]
%    ds_stim     The downsampled, preprocessed stimulus   [SxTxF]
%    ds_resp     The downsampled, preprocessed response   [SxTxRxF]
%    ds_time     The downsampled time index               [1xT]
%    ds_psth     The downsampled, averaged response       [SxTxF]
%
% In the above, dimensions are indicated with
%      S = sound stimulus index #
%      R = repetition index #
%      N = Time index at 100KHz sampling
%      T = Time index in downsampled frequency
%      F = Preprocessing index #

log_dbg('load_stim_resps() was called');

rasterfreq = 100000;   % TODO: Make user editable
includeprestim = 1;    % TODO: Make user editable

% Merge the training and test set names, which may overlap
files_to_load = unique({training_set{:}, test_set{:}});

% Create the 'dat' cell array and its entries
len = length(files_to_load);
dat = cell(1, len);
for i = 1:len;
    f = files_to_load{i};
    
    % Find the right index of cfd corresponding to the file, call it idx
    idx = 0;
    for j = 1:length(cfd)
        if isequal(cfd(j).stimfile, f)
            idx = j;
        end
    end
    if idx == 0 log_err('Not found in cfd: %s', f); end
    
    % Load the raw_stim part of the data structure
    stimfile = [cfd(idx).stimpath cfd(idx).stimfile];
    log_msg('Loading stimulus: %s', stimfile);
    stim = loadstimfrombaphy(stimfile, [], [], 'wav', rasterfreq, 1, 0,...
                             includeprestim);
    [d1 d2 d3] = size(stim);
    if d1 ~= 1
         log_dbg('Stimulus matrix was: [%d %d %d]', d1, d2, d3);
         log_err('Stimulus size was not [1xNxS]!?');
    end
    dat.(f).raw_stim = permute(stim, [3 2 1]); 
    
    % Create the raw_time index
    dat.(f).raw_time = (1/rasterfreq).*[1:d2]';
    
    % Load the raw_resp part of the data structure
    respfile = [cfd(idx).path, cfd(idx).respfile];
    log_msg('Loading response: %s', respfile);
    options = [];
    options.includeprestim = includeprestim;
    options.unit = cfd(idx).unit;
    options.channel  = cfd(idx).channum;
    options.rasterfs = rasterfreq; 
    [resp, tags] = loadspikeraster(respfile, options);   
    dat.(f).raw_resp = permute(resp, [3 1 2]);    
    dat.(f).raw_respavg = squeeze(sum(dat.(f).raw_resp, 3));
    
    % Check raw_stim, raw_resp, and raw_time signal sizes match.
    [s1 s2]    = size(dat.(f).raw_stim);
    [r1 r2 r3] = size(dat.(f).raw_resp);
    [a1 a2]    = size(dat.(f).raw_respavg);
    if ~(isequal(s1, r1, a1) & isequal(s2, r2, a2))
        log_dbg('Stim [%d %d], Resp: [%d %d %d] Respavg=[%d %d]',...
                s1,s2,r1,r2,r3, a1,a2);
        log_err('Stimulus, Response, and Average matrices size mismatch.');
    end
end
log_dbg('Done loading stimulus and response files.');

% ------------------------------------------------------------------------
function fns = scan_directory_for_functions(scanpath)
fprintf('Scanning dir for functions: %s\n', scanpath);
files = dir(fullfile(scanpath, '*.m'));
fns = [];
for i = 1:length(files)
    name = files(i).name;
    fprintf('\tFound ''%s''\n', name);
    s = [];                        % The struct to save
    s.fn_name = name(1:end-2);     % File name minus the '.m' at end
    fn = str2func(name(1:end-2));  % Make an executable function
    s.fn = fn;                     % Save function handle
    s.params = fn();               % Default param struct given by no args
    fns.(s.fn_name) = s;           % Index under its function name
end

% ------------------------------------------------------------------------
function s = repl_write(obj)
% Prints obj in a readable manner. 
if ismatrix(obj) & isnumeric(obj) & any(size(obj) ~= 1)  % Matrices
    s = strcat('[', num2str(obj), ']');
    s = regexprep(s, '\n', '; ');
elseif isstr(obj)       % Single strings
    s = obj; 
elseif isnumeric(obj)   % Single numbers
    s = num2str(obj);
elseif isa(obj, 'function_handle')
    s = ['@' func2str(obj)];
else
    log_err('Not sure how to print: %s', obj);
end

%------------------------------------------------------------------------
function s = extract_field_val_pairs(mytable, fieldname_col, value_col)
% Return a new struct extracted from two columns
d = get(mytable, 'Data');
[r, c] = size(d);
if (c < fieldname_col | c < value_col | fieldname_col < 1 |  value_col < 1)
    err('Column index number is outside the data table''s range.');
end
s = {};
for i = 1:r
    s.(d{i,fieldname_col}) = eval(d{i,value_col});
end

%------------------------------------------------------------------------
function s = extract_checked_fields(mytable, checkbox_col, fieldname_col)
% Return a cell array of fields with checked boxes next to them.
d = get(mytable, 'Data');
[r, c] = size(d);
j = 1;
s = {};
for i = 1:r
    if d{i,checkbox_col}
        s{j} = d{i,fieldname_col};
        j = j+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA SELECTION GUI

function query_db_button_Callback(hObject, eventdata, handles)
log_dbg('query_db_button pressed.');

% GS: Clear various parts of the global struct
global GS;

% Invalidate data_selection_table, and other parts of the GUI
set(handles.data_selection_table, 'Data', {}); drawnow;
axes(handles.stim_view_axes); cla;
axes(handles.resp_view_axes); cla;
axes(handles.preproc_view_axes); cla;
axes(handles.model_view_axes); cla;
axes(handles.optplot1); cla;
axes(handles.optplot2); cla;
axes(handles.optplot3); cla;
set(handles.selected_stimfile_popup, 'String', '');   
set(handles.selected_stim_idx_popup, 'String', '');  

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
function update_raw_stim_plot(handles)
log_dbg('Updating raw_stim_plot');
global GS;
 
% Only update if all fields are defined
if ~(isfield(GS, 'raw_stim_plot_type') & ...
     isfield(GS, 'selected_stim_idx') & ...
     isfield(GS, 'selected_stimfile'))  
     log_dbg('Ignoring raw stim plot update since not all fields ready.');
   return  
end

plottype = GS.raw_stim_plot_type;
stim_idx = GS.selected_stim_idx;
stimfile = GS.selected_stimfile;
obj = GS.dat.(stimfile);
   
axes(handles.stim_view_axes);cla;
switch plottype
    case 'Time Series View'
        plot(obj.raw_time, obj.raw_stim(stim_idx,:), 'k-');
        axis tight;
    case 'Spectrogram View'
        % From 500Hz, 12 bins per octave, 4048 sample window w/half overlap
        logfsgram(obj.raw_stim(stim_idx,:)', 4048, 100000, [], [], 500, 12); 
        % TODO: Remove 4048, 100000 here and use global
        caxis([-20,40]);  % TODO: use a 'smarter' caxis here
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
preprocessing_view_popup_Callback(handles.preprocessing_view_popup, [], handles);
% TODO: update_model_view_plot(handles);

%------------------------------------------------------------------------
function selected_stim_idx_popup_CreateFcn(hObject, eventdata, handles)
function selected_stim_idx_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(handles.selected_stim_idx_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stim_idx_popup, 'Value')});
GS.selected_stim_idx = stim_idx;
update_raw_stim_plot(handles);
update_raw_resp_plot(handles);
preprocessing_view_popup_Callback(handles.preprocessing_view_popup, [], handles);
% TODO: update_model_view_plot(handles);

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
function update_data_table(mytable, mystruct, myfields)
l = length(myfields);
c = cell(l,3);
for i = 1:l
    if ~isfield(mystruct, myfields{i})
        log_err('Could not find field: %s', myfields{i});
    end
    c{i,1} = false;
    c{i,2} = myfields{i};
    c{i,3} = repl_write(mystruct.(myfields{i})); % Ensure data becomes a str
end
set(mytable, 'Data', c);
drawnow;

%------------------------------------------------------------------------
function preproc_popup_CreateFcn(hObject, eventdata, handles)
function preproc_popup_Callback(hObject, eventdata, handles)
global GS; 
c = cellstr(get(hObject,'String'));
pretty_name = c{get(hObject,'Value')};

% Find the selected preproc filter
% TODO: Replace me with find() idiom instead of a fucking for loop
f = [];
fnames = fieldnames(GS.preprocs);
for i = 1:length(fnames)  
    if isequal(pretty_name, GS.preprocs.(fnames{i}).params.pretty_name)
        f = GS.preprocs.(fnames{i});
    end
end

% If f is not found throw an error
if ~isequal(f, [])
    GS.selected_preproc_name = f.fn_name;
else
	log_err('Somehow, the selected preproc name was not found!?');
end
% --------------------

% TODO: Replace with two function calls that set values and checkboxes independently?
pp = GS.preprocs.(GS.selected_preproc_name).params;
update_data_table(handles.preproc_data_table, pp, pp.editable_fields);
              
%------------------------------------------------------------------------
function apply_preproc_button_Callback(hObject, eventdata, handles)
global GS;
log_inf('Preprocessing...');

% TODO: Move this core computation somewhere else!
f = fieldnames(GS.dat); 
spp = GS.preprocs.(GS.selected_preproc_name).params;
for i = 1:length(f)
    % Feed data into that new filter
    GS.dat.(f{i}).pp_stim = spp.preproc_fn(spp, GS.dat.(f{i}).raw_stim); 
end

% TODO: Check that the length of the preprocessed vector is the SAME size
% as the original raw stimulus.

% Call the callbacks which initialize thingscallback
set(handles.preproc_index_popup, 'String', 'NONE');
preproc_index_popup_Callback(handles.preproc_index_popup, [], handles);

log_inf('Done preprocessing.');
preprocessing_view_popup_Callback(handles.preprocessing_view_popup, [], handles);

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

spp = GS.preprocs.(GS.selected_preproc_name);

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
            axis tight;
        end
     case 'Filtered Spectrogram'
        set(handles.preproc_index_popup, 'Enable', 'On');
        set(handles.view_preproc_index_label, 'Enable', 'On');
        logfsgram(dat.pp_stim(GS.selected_stim_idx,:,GS.selected_preproc_idx)', 4048, 100000, [], [], 500, 12); 
        caxis([-20,40]);
        drawnow;
end
hold off;

%------------------------------------------------------------------------
function preprocessing_view_popup_CreateFcn(hObject, eventdata, handles)
function preprocessing_view_popup_Callback(hObject, eventdata, handles)
global GS;
nvals = cellstr(get(hObject, 'String'));
plottype = nvals{get(hObject, 'Value')};
GS.preproc_view_plot_type = plottype;
update_preproc_view_plot(handles);

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

%------------------------------------------------------------------------
function load_preproc_params_button_Callback(hObject, eventdata, handles)
function save_preproc_params_button_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------------
function preproc_data_table_CellEditCallback(hObject, eventdata, handles)
global GS;
log_dbg('preproc_data_table modified. Click ''preprocess!'' to refresh.');
% Whenever the data table is edited, do three things:
% 1. Pull out the present values and update the params struct first
s = extract_field_val_pairs(hObject, 2, 3);
fns = fieldnames(s);
pp = GS.preprocs.(GS.selected_preproc_name);
for i = 1:length(fns);
    pp.params.(fns{i}) = s.(fns{i});
end
GS.preprocs.(GS.selected_preproc_name).params = pp.fn(pp.params); 

% 2. Update which parameters are desired to be fit with the optimization
GS.preprocs.(GS.selected_preproc_name).fittable_params = ...
     extract_checked_fields(hObject, 1, 2);

% 3. Invalidate by setting preprocessed data to [] and clearing plot
GS.dat.(GS.selected_stimfile).pp_stim = [];
axes(handles.preproc_view_axes); cla

% --- Executes when entered data in editable cell(s) in preproc_data_table.
function auto_preproc_data_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to preproc_data_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOWNSAMPLING WIDGETS

function downsample_stimresp()
global PF_STIM RESPAVG SAMPFREQ DS_FREQ DS_STIM DS_TIME DS_RESPAVG FIRBINSIZE;
DS_FREQ = 1000 / FIRBINSIZE;  
disp('Downsampling the PF_STIM and RESP signals...');
DS_RESPAVG = conv_fn(RESPAVG, 2, @sum, SAMPFREQ/DS_FREQ, 0);
[d,l] = size(DS_RESPAVG);
DS_TIME = [0:1/DS_FREQ:l/DS_FREQ-1/DS_FREQ];
DS_STIM = conv_fn(PF_STIM, 2, @mean, SAMPFREQ/DS_FREQ, 0);

% TODO: 
%    case 'Downsampled Stimulus'
%        % If the downsampled data exists, plot it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CLASS WIDGETS
function update_fir_model(handles)
global PF_COEFS FIRHIST FIRBINSIZE FIRCOEFS;
FIRHIST = str2num(get(handles.fir_history,'String'));
FIRBINSIZE = str2num(get(handles.bin_size,'String'));
d = length(PF_COEFS); % Since this is a cell array, this returns just nfilts
l = ceil(FIRHIST/FIRBINSIZE);
FIRCOEFS = zeros(d,l);

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
    DS_PREDS = cat(3, DS_PREDS, sqrt(abs(filter(FIRCOEFS(filt_idx,:), ...
        1, squeeze(DS_STIM(:,:,filt_idx)), [],2))));
end

DS_PRED = squeeze(sum(DS_PREDS, 3)); 

function plot_model(handles)
global FIRCOEFS DS_TIME DS_STIM DS_RESPAVG DS_PREDS DS_PRED FIRBINSIZE;

nvals = cellstr(get(handles.model_view_popup, 'String'));
plottype = nvals{get(handles.model_view_popup, 'Value')};

nvals = cellstr(get(handles.selected_stim_idx_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stim_idx_popup, 'Value')});

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

nvals = cellstr(get(handles.selected_stim_idx_popup, 'String'));
stim_idx = str2num(nvals{get(handles.selected_stim_idx_popup, 'Value')});

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


function downsampling_popup_CreateFcn(hObject, eventdata, handles)
function downsampling_popup_Callback(hObject, eventdata, handles)
function load_downsampling_params_button_Callback(hObject, eventdata, handles)
function save_downsampling_params_button_Callback(hObject, eventdata, handles)

function load_model_params_button_Callback(hObject, eventdata, handles)
function save_model_params_button_Callback(hObject, eventdata, handles)
function predict_response_button_Callback(hObject, eventdata, handles)

function load_stochasticity_params_button_Callback(hObject, eventdata, handles)
function save_stochasticity_params_button_Callback(hObject, eventdata, handles)
function check_stochasticity_button_Callback(hObject, eventdata, handles)

function pushbutton22_Callback(hObject, eventdata, handles)


% --- Executes on selection change in popupmenu21.
function popupmenu21_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu21 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu21


% --- Executes during object creation, after setting all properties.
function popupmenu21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
