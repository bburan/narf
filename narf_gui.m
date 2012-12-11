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

% Last Modified by GUIDE v2.5 26-Nov-2012 13:41:07

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
narf_set_path;
global NARF_PATH NARF_MODULES_PATH NARF_OPTIMIZATION_PATH;

%%%%%%%%%%%%%% GLOBAL STUFF ABOVE %%%%%%%%%%%%%%%%%%

% Build a list of the files in the 'modules/' directory
mods = scan_directory_for_modules(NARF_MODULES_PATH);

% Add the model pane
global STACK XXX;
STACK = {};
XXX = {};
GS.refresh_gui = narf_modelpane(handles.model_structure_panel, mods);
    
% % Invalidate all data tables
set(handles.data_selection_table, 'Data', {});
set(handles.sampling_data_table, 'Data', {});
set(handles.perf_metric_data_table, 'Data', {});
set(handles.term_cond_data_table, 'Data', {});
drawnow;

% Set the default cellid to be something relevant
set(handles.cellid_text, 'String', 'por022b-a1');

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
% parms = cell(1, len);
% perfs = cell(1, len);

% Load parms, perfs, if you need them
% [parms{i}, perfs{i}] = dbReadData(cfd(i).rawid);

% Select set with the most repetitions as test set.
for i = 1:len;
    if (isfield(cfd(i), 'repcount') && cfd(i).repcount > test_set_reps)
       test_set_reps = dbget('gDataRaw', cfd(i).rawid, 'reps');
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

% TODO: Remove me later, this is just for making things happen faster while
% I test the same functions over and over again
training_set=training_set(1);

% TODO: Consider saving the parm/perf info to a global for later use?
% TODO: Consider moving parm/perf outside this fn?
log_dbg('Training sets selected: %s', ...
    strtrim(sprintf('%s ', training_set{:})));
log_dbg('Test set selected: %s', char(test_set{:}));
 
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
global GS STACK XXX;
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

% Reset the values of the XXX structure
XXX = {};
XXX{1} = [];
XXX{1}.cellid = GS.cellid;
XXX{1}.training_set = GS.training_set;
XXX{1}.test_set = GS.test_set;

% Request that the model recompute values, if they exist
if isfield(GS, 'refresh_gui')
    GS.refresh_gui();
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOWNSAMPLING WIDGETS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CLASS WIDGETS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STOCHASTICITY AND ITS PLOTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMIZATION WIDGETS

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
global STACK;
% Get the number of iterations
n_iters = 10;

% Get the starting score of the current filter
x_0 = pack_fittables(STACK);

if isempty(x_0)
    log_msg('No parameters were selected to be fit.');
    return;
end

stepsize = 1.0;
[x_bst, s_bst] = boosting(x_0', @correlation_of_downsampled_signals, @(n,x,s)(n > n_iters), stepsize);

unpack_fittables(x_bst);

% Finally, update the GUI's data tables and plots 
update_tables_and_plots();

function initialize_optplot_menu(handle)
menuopts = {'Pred/Resp Scatter', 'Raw ISI', 'Intensity-Scaled ISI', ...
        'KS Plot'};
set(handle, 'String', char(menuopts));

function setoptplot(handles, plothandle, plottype)
global GS;

% stim_idx = GS.selected_stim_idx;
% dat = GS.dat.(GS.selected_stimfile);
% 
% axes(plothandle); cla;
% hold on;
% switch plottype
%     case 'Pred/Resp Scatter'
%         plot(dat.ds_respavg(stim_idx,:), dat.ds_pred(stim_idx,:), 'k.');
%         axis tight;
%     case 'Raw ISI'
%         %plot(DS_PRED, DS_RESPAVG, 'k.');
%         %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
%         axis tight;        
%     case 'Time-Scaled ISI'
%         %plot(DS_PRED, DS_RESPAVG, 'k.');
%         %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
%         axis tight;
%     case 'KS Plot'
%         %plot(DS_PRED, DS_RESPAVG, 'k.');
%         %setAxisLabelCallback(gca, @(t) (FIRBINSIZE*(t-1)), 'X');
%         axis tight;
% end
% hold off;

function optplot1popup_CreateFcn(hObject, eventdata, handles)
initialize_optplot_menu(hObject);

function optplot1popup_Callback(hObject, eventdata, handles)
nvals = cellstr(get(handles.optplot1popup, 'String'));
plottype = nvals{get(handles.optplot1popup, 'Value')};
setoptplot(handles, handles.optplot1, plottype);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_stochast_button_Callback(hObject, eventdata, handles)

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
