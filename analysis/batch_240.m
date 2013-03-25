function batch_240()

% MODEL BLOCK GROUPS THAT WORK WELL (as of 2013-03-06)
% 
% mm{1}: LOADER GROUP        Gets the data from Baphy
%
% mm{2}: COMPRESSOR GROUP    Accounts for logarithmic sensitivity to volume
%        log2b: log2 with baseline so output is positive semidefinite.
%        log2:  log2 without baseline (so may go negative)
%
% mm{3}: FILTER BLOCK        Accounts for time latency
%        firn: Normalization on inputs, normalization on output
%        depn: 13 channel depression filter bank with normalization.
%
% mm{4}: NONLINEARITY GROUP  Accounts for input-output nonlinearity
%        npnl: Binned means. Good, fast, but extrapolates poorly.
%        npfnl: Filter-based npnl five gaussians wide.
%        npfnl4: Filter-based npnl five gaussians wide.
%        npfnl3: Filter-based npnl three gaussians wide. 
%        senl: Sparse gaussian mixture centered at 'representative' points.
%        senl25: Same, but with slightly wider gaussians (0.25 width)
%        senl3: Same, but with slightly wider gaussians (0.3 width)
%
% mm{5}: FITTING GROUP       Determines how the model is fit
%        boost: A simple boosting implementation.
%        sb: Sparse bayesian boosting with 'sparse' steps
%        fmin: Basic line search built in to matlab
%        lsq: Basic lsqcurvefit() built into matlab.
%        fminlsq: First fmin, then lsq
%
% mm{6}: OBJECTIVE FUNCTION  Measures performance, applies penalties, 
%        mse  : Performance metric is mean squared error, with no smoothness or sparseness penalty.
%        mses1: Same as MSE, but a sparseness penalty of 10^-1. Strong sparseness penalty
%        mses2: Same as MSE, but a sparseness penalty of 10^-2.
%        mses3: Same as MSE, but a sparseness penalty of 10^-3. 
%        mses4: Same as MSE, but a sparseness penalty of 10^-4. 
%        mses5: Same as MSE, but a sparseness penalty of 10^-5. 
%        mses6: Same as MSE, but a sparseness penalty of 10^-6. Weak sparseness penalty  

% CAVEATS:
%   'sb' fitter and 'senl' nonlinearities crash about once every few
%   hundred models because of numerical instabilities. They work really
%   well when they work, though.
batch = 240;
analysis_prefix = 'a240';
baphy_set_path;
narf_set_path; 

global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();   

cells = request_celldb_batch(batch);

mm = {};
mm{1} = module_groups('env100');
mm{2} = module_groups('log2b'); 
mm{3} = module_groups('firn');  
mm{4} = module_groups('npfnl');   % 'npnl', 'senl', 'senl3',
mm{5} = module_groups('sb', 'fminlsq', 'boost', 'fmin');
mm{6} = module_groups('mse', 'mses2', 'mses3', 'mses4', 'mses5', 'mses6'); 

modulekeys = module_block_combos(mm);

force = false;

% Enqueue every model
for ii = 1:length(cells)
    for jj = 1:length(modulekeys)
        %fprintf('Fitting model [%d/%d]\n', jj, length(modulekeys)); 
        %fit_single_model(modulekeys{jj}, batch, cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
        % TODO: Enqueue into job system instead of doing fit_single_model
        % here to allow work to be distributed everywhere.
        
        force = true;
        
        enqueue_single_model(modulekeys{jj},  batch, cells{ii}.cellid, ...
            cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, force);
        
    end
end

% Generate "top 10" plots for each cellid
% for ii = 1:length(cells)
%     % Plot the top 10 models
%     models = db_get_models(batch, cells{ii}.cellid);
%     compare_models(cellstr(char(models(1:min(length(models), 10)).modelpath)));
%     % TODO: Scatter plots
% end

% TODO: Generate heat map plots (Buggy right now due to caching problem)


% Write a function which makes an SQL query
% Look at db_get_models() for an example
%

