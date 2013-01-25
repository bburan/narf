function results = test_likely_candidates(cellid, training_set, test_set)
% Exhaustively test the most likely model structures candidates for fitting
% This is the "throws the entire kitchen at the problem, including the 
% sink" approach to determining which model describes the activity of a
% particular neuron.
%
% Arguments are self explanatory.
% 
% Returns a long cell array of structs.
%
% EXAMPLE USE:
%  test_likely_candidates('por028d-b1', {'por029d03_p_SPN'}, {'por028d04_p_SPN'});

global STACK XXX MODULES NARF_PATH NARF_SAVED_MODELS_PATH ...
    NARF_MODULES_PATH NARF_SAVED_ANALYSIS_PATH;

if isempty(MODULES)
    MODULES = scan_directory_for_modules(NARF_MODULES_PATH);
end

mm = {}; 

% GROUP 1: LOADING A DOWNSAMPLED, POSITIVE SEMIDEFINITE SIGNAL
mm{1} = [];
% mm{1}.env200hz = {MODULES.load_stim_resps_from_baphy.mdl(...
%                                     struct('raw_resp_fs', 200, ...
%                                            'raw_stim_fs', 200,...
%                                            'stimulus_format', 'envelope'))};
mm{1}.env100hz = {MODULES.load_stim_resps_from_baphy.mdl(...
                                    struct('raw_resp_fs', 100, ...
                                           'raw_stim_fs', 100,...
                                           'stimulus_format', 'envelope'))};                                    
% mm{1}.wav100khz_elliptic_ds200hz = ...
%                  {MODULES.load_stim_resps_from_baphy.mdl(...
%                                     struct('raw_resp_fs', 200, ...
%                                            'raw_stim_fs', 100000)), ...
%                   MODULES.elliptic_bandpass_filter_bank, ...
%                   MODULES.downsample_with_fn.mdl(...
%                                     struct('downsampled_freq', 200, ...
%                                            'postconv_fn', @(x) x))};
% mm{1}.wav100khz_elliptic_ds100hz = ...
%                  {MODULES.load_stim_resps_from_baphy.mdl(...
%                                     struct('raw_resp_fs', 200, ...
%                                            'raw_stim_fs', 100000)), ...
%                   MODULES.elliptic_bandpass_filter_bank, ...
%                   MODULES.downsample_with_fn.mdl(...
%                                     struct('downsampled_freq', 200, ...
%                                            'conv_fn', @mean, ...
%                                            'postconv_fn', @(x) x))};

% TODO: I needed the abs() around the sqrt() because for some reason
% load_stim_from_baphy is giving negative zero values for the envelope,
% which then turn complex and crash the sigmoidal nonlinearity?

% GROUP 2: COMPRESSION OF INPUT INTENSITY
mm{2} = [];
mm{2}.nocomp = {MODULES.passthru};
% The nth-root compressors do not work as well as log compressors, but if
% you want to try them feel free.
% mm{2}.root15 = {MODULES.nonlinearity.mdl(struct('phi', [], ...
%                                                'nlfn', @(phi, z) abs(z.^(1/1.5))))};
% mm{2}.root2 = {MODULES.nonlinearity.mdl(struct('phi', [], ...
%                                                'nlfn', @(phi, z) abs(sqrt(z))))};
% mm{2}.root25 = {MODULES.nonlinearity.mdl(struct('phi', [], ...
%                                                'nlfn', @(phi, z) abs(z.^(1/2.5))))};
% mm{2}.root3 = {MODULES.nonlinearity.mdl(struct('phi', [], ...
%                                                 'nlfn', @(phi, z) abs(z.^(1/3))))};
% mm{2}.root4 = {MODULES.nonlinearity.mdl(struct('phi', [], ...
%                                                 'nlfn', @(phi, z) z.^(1/4)))};
% mm{2}.root5 = {MODULES.nonlinearity.mdl(struct('phi', [], ...
%                                                'nlfn', @(phi, z) z.^(1/5)))};
% mm{2}.rootfit = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
%                                                  'phi', [1/2.5], ...
%                                                  'nlfn', @(phi, z) abs(z.^(phi(1)))))};

% Over a population of 33 cells, some model in the range log1 to log4 was
% the best, so I don't think we need to search outside that. 
mm{2}.log1  = {MODULES.nonlinearity.mdl(struct('phi', [], ...
                                               'nlfn', @(phi, z) log(z + 10^-1)))};
mm{2}.log2  = {MODULES.nonlinearity.mdl(struct('phi', [], ...
                                               'nlfn', @(phi, z) log(z + 10^-2)))};
mm{2}.log3  = {MODULES.nonlinearity.mdl(struct('phi', [], ...
                                               'nlfn', @(phi, z) log(z + 10^-3)))};
mm{2}.log4  = {MODULES.nonlinearity.mdl(struct('phi', [], ...
                                               'nlfn', @(phi, z) log(z + 10^-4)))};
%mm{2}.log5  = {MODULES.nonlinearity.mdl(struct('phi', [], ...
%                                               'nlfn', @(phi, z) log(z + 10^-5)))};

% mm{2}.logfit  = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
%                                                  'phi', [-13], ...
%                                                  'nlfn', @(phi, z) log(z + exp(phi(1)))))};

%mm{2}.volterra = {MODULES.concat_second_order_terms}; 
%mm{2}.depress  = {MODULES.depression_filter_bank}; 

% GROUP 3: THE FIR FILTER 
mm{3} = [];
% mm{3}.fir =      {MODULES.normalize_channels, ...
%                   MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
%                                                 'fit_fields', {{'coefs'}}))};
mm{3}.firbase =  {MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs', 'baseline'}}))};

                                            
                                            %mm{3}.depfir =  {MODULES.depression_filter_bank, ...
%                 MODULES.normalize_channels, ...
%                 MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
%                                               'fit_fields', {{'coefs', 'baseline'}}))};

% 
% mm{3}.inhibexcit = {MODULES.normalize_channels, ...
%                     MODULES.fir_filter.mdl(struct('fit_fields', {{'coefs'}}, ...
%                                                   'num_coefs', 12, ...
%                                                   'output', 'inhib')), ...
%                     MODULES.fir_filter.mdl(struct('fit_fields', {{'coefs'}}, ...
%                                                   'num_coefs', 12, ...
%                                                   'output', 'excit')), ...                                                  
%                     MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
%                                                     'input_stim', 'inhib', ...
%                                                     'output', 'inhib', ...
%                                                     'phi', [0], ...
%                                                     'nlfn', @(phi, z) - zero_below_thresh(phi, z))), ...
%                     MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
%                                                     'input_stim', 'excit', ...
%                                                     'output', 'excit', ...
%                                                     'phi', [0], ...
%                                                     'nlfn', @zero_below_thresh)), ...
%                     MODULES.sum_fields.mdl(struct('inputs', {{'inhib', 'excit'}}, ...
%                                                   'output', 'stim'))};

% GROUP 4: POST-FILTER NONLINEARITIES
mm{4} = [];
mm{4}.nonl = {MODULES.passthru};
mm{4}.sig  = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                             'phi', [0 1 1 0], ...
                                             'nlfn', @sigmoidal))};
mm{4}.npnl = {MODULES.nonparm_nonlinearity};

% mm{4}.exp  = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
%                                               'phi', [1 1], ...
%                                               'nlfn', @exponential))};
% mm{4}.step = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
%                                               'phi', [0], ...
%                                               'nlfn', @zero_below_thresh))};
% mm{4}.poly = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
%                                               'phi', [1 1 1 1], ...
%                                               'nlfn', @polyval))};       

% GROUP 5: PERFORMANCE METRICS AND OPTIMIZATION TECHNIQUES
% TODO: Right now the names of these are critical because there is no way
% of hinting to the optimization routines how they should be tested.

% In general, fminlsq seems to work the best.
mm{5} = [];
mm{5}.fmin = {MODULES.correlation, ...
              MODULES.mean_squared_error.mdl(struct('output', 'score'))};   
mm{5}.lsq  = {MODULES.correlation, ...
              MODULES.mean_squared_error.mdl(struct('output', 'score'))};  
mm{5}.fminlsq  = {MODULES.correlation, ...
              MODULES.mean_squared_error.mdl(struct('output', 'score'))};
mm{5}.smooth = {MODULES.correlation, ...
               MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                              'smoothness_weight', 10^-6))};
mm{5}.jack = {MODULES.correlation, ...
               MODULES.mean_squared_error.mdl(struct('output', 'score'))};
mm{5}.twostep = {MODULES.correlation, ...
                 MODULES.mean_squared_error.mdl(struct('output', 'score'))};
mm{5}.threestep = {MODULES.correlation, ...
                 MODULES.mean_squared_error.mdl(struct('output', 'score'))};
% ------------------------------------------------------------------------
% BUILD THE MODELS

opts = cellfun(@fieldnames, mm, 'UniformOutput', false);
N_opts = cellfun(@(m) length(fieldnames(m)), mm);
N_models = prod(N_opts);
fprintf('Number of models to be tested: %d\n', N_models); 

results = cell(N_models, 1);
mkdir([NARF_SAVED_MODELS_PATH filesep cellid]);
analysis_file = [NARF_SAVED_ANALYSIS_PATH filesep cellid '_results.mat'];

% Load existing analysis file, so that incomplete analyses may be continued
if exist(analysis_file, 'file') == 2
    fprintf('Loading existing analysis file.\n');
    nnn = load(analysis_file, 'results');
    results = nnn.results;
end

for ii = 1:N_models,
    STACK = {};
    XXX = {};
    XXX{1}.cellid = cellid;
    XXX{1}.training_set = training_set;
    XXX{1}.test_set = test_set;
    
    tic;
    
    % Behold matlab's ugliness for destructuring binds!
    [i1 i2 i3 i4 i5 i6 i7 i8 i9] = ind2sub(N_opts, ii);
    indexes = [i1 i2 i3 i4 i5 i6 i7 i8 i9];
    
    opt_names = arrayfun(@(gi, ind) opts{gi}{ind}, ...
        1:length(opts), indexes(1:length(opts)), ...
        'UniformOutput', false);
    
    blocks = cellfun(@getfield, mm, opt_names, 'UniformOutput', false);
    
    tmp = cellfun(@(n) sprintf('%s_', n), opt_names, 'UniformOutput', false);
    modelname = strcat(tmp{:});
    modelname = modelname(1:end-1); % Remove last underscore character
    modelfile = [NARF_SAVED_MODELS_PATH filesep cellid filesep cellid '_' modelname '_' strcat(training_set{:}) '.mat'];
    fprintf('MODEL [%d/%d]: %s\n', ii, N_models, modelname);
    
    % If the model savefile exists, we assume we don't need to fit 
    if exist(modelfile, 'file') == 2
        fprintf('Skipping because model file exists.\n');
        continue;
    end
    
    % Append each block one at a time
    for jj = 1:length(blocks),
        fprintf('Auto-initalizing module [%d/%d]\n', jj, length(blocks));
        for kk = 1:length(blocks{jj})
            append_module(blocks{jj}{kk}); % Calls auto_init for each
        end
    end
    
    % TODO: 
    % The correct way to fit things is to realize that each fitter will
    % work better in different circumstances, and there is no perfect
    % fitting routine which works for everything. 
    % Therefore, we simply try several fit routines, and pick the
    % best-performing parameters found among all fit routines. 

    fprintf('Fitting all parameters\n');
    if strcmp(opt_names{end}, 'fmin'),
        exit_code = fit_objective('score');
    elseif strcmp(opt_names{end}, 'fminlsq')
        fit_objective('score');
        exit_code = fit_with_lsqcurvefit();
    elseif strcmp(opt_names{end}, 'lsq')
        exit_code = fit_with_lsqcurvefit();
    elseif strcmp(opt_names{end}, 'smooth')
        exit_code = fit_objective('score');
    elseif strcmp(opt_names{end}, 'jack')
        exit_code = fit_with_jacklsq();
    elseif strcmp(opt_names{end}, 'twostep')
        % Fit FIR by itself, then NL by itself
        indexes = find_module_indexes(STACK, 'nonlinearity');
        [firmod, firmodidx] = find_module(STACK, 'fir_filter');
        if ~isempty(indexes) && isfield(STACK{indexes(end)}, 'fit_fields')
            firfields = firmod.fit_fields;
            nlfields = STACK{indexes(end)}.fit_fields;
            % Step 1: Just FIR
            STACK{indexes(end)}.fit_fields = {};
            fit_objective('score');
            fit_with_lsqcurvefit();
            % Step 2: Just NL
            STACK{firmodidx}.fit_fields = {};
            STACK{indexes(end)}.fit_fields = nlfields;
        end
        fit_objective('score');
        exit_code = fit_with_lsqcurvefit();
    elseif strcmp(opt_names{end}, 'threestep')
        % Fit FIR by itself, then NL by itself
        indexes = find_module_indexes(STACK, 'nonlinearity');
        [firmod, firmodidx] = find_module(STACK, 'fir_filter');
        if ~isempty(indexes) && isfield(STACK{indexes(end)}, 'fit_fields')
            firfields = firmod.fit_fields;
            nlfields = STACK{indexes(end)}.fit_fields;
            % Step 1: Just FIR
            STACK{indexes(end)}.fit_fields = {};
            fit_objective('score');
            fit_with_lsqcurvefit();
            % Step 2: Just NL
            STACK{firmodidx}.fit_fields = {};
            STACK{indexes(end)}.fit_fields = nlfields;
            fit_objective('score');
            fit_with_lsqcurvefit();
            % Step 3: BOTH
            STACK{firmodidx}.fit_fields = firfields;
        end
        fit_objective('score');
        exit_code = fit_with_lsqcurvefit();
    end
    
    % TODO: Move this to 'make_analysis_cache'
    results{ii}.score_train_corr = XXX{end}.score_train_corr;
    results{ii}.score_test_corr = XXX{end}.score_test_corr;
    results{ii}.score_train_mse = XXX{end}.score_train_mse;
    results{ii}.score_test_mse = XXX{end}.score_test_mse;
    firmod = find_module(STACK, 'fir_filter');
    results{ii}.fir_coefs = firmod.coefs;    
    indexes = find_module_indexes(STACK, 'nonlinearity');
    if ~isempty(indexes)
        results{ii}.nl_phi = STACK{indexes(end)}.phi;
    end
    results{ii}.fit_time = toc;
    results{ii}.n_free_params = length(pack_fittables(STACK));
    results{ii}.exit_code = exit_code;    
    results{ii}.cellid = XXX{1}.cellid;
    results{ii}.training_set = XXX{1}.training_set;
    results{ii}.test_set = XXX{1}.test_set;
    results{ii}.optimized_on = STACK{end}.name;
    results{ii}.filename     = modelfile;
    results{ii}.modelname    = modelname;
    
    % TODO: Remove this when I'm done debugging
    XXX{1}.results_cache = results{ii};
    
    save_model_stack(modelfile, STACK, XXX);
    
    % Save results incrementally to analysis file
    save(analysis_file, 'results');
    
end
