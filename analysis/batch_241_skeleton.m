function batch_241_skeleton()

% LIMITATIONS:
%   - Single machine, single matlab instance only 
%   - Easy to generate many models, but hard to aggregate data from all of them
% 
% Both of these will be fixed ASAP. 

batch = 241;
baphy_set_path;
narf_set_path; 

global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();
    
analysis_prefix = 'a241';
    
cells = request_celldb_batch(batch);

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
%   well when the work, though.

% Recommended defaults for minimal search. As desired, add the above
mm = {};
mm{1} = module_groups('env100');
mm{2} = module_groups('log2b'); 
mm{3} = module_groups('firn');  
mm{4} = module_groups('npnl');  % senl and npfnl4 usually beat npnl
% TODO: NPNLX module group, which uses separate NPNLs for each respfile
mm{5} = module_groups('sb', 'fminlsq', 'boost');  % The 3 best fitters
mm{6} = module_groups('mse', 'mses1', 'mses5'); % Zero, heavy, and light sparsity
    
[~, modelnames] = module_combinations(mm);

% Generate every possible model combination for each cellid in this batch
for ii = 1:length(cells)
    fit_models(mm, cells{ii}.cellid, cells{ii}.training_set, cells{ii}.test_set);
end

% Generate "top 10" plots and scatter plots for each cellid
for ii = 1:length(cells)
    % ----------- TO BE REPLACED WITH MYSQL CODE --------
    % The following will hopefully be irrelevant soon because its SUPER DUPER DOG SLOW
    summarize_cellid(cells{ii}.cellid, true);
    sf = [NARF_SAVED_ANALYSIS_PATH filesep cells{ii}.cellid '_summary.mat'];
    summaries = load_summaries({sf});
    summaries = only_named_summaries(summaries, modelnames);
    plot_cellid_summary(cells{ii}.cellid, summaries, true, analysis_prefix);
    % -------------------- END REPLACEMENT -------------- 
end

% TODO: Generate heat map plots (Buggy right now due to caching problem)