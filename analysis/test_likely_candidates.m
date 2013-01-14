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
%  test_likely_candidates('por022a-a1', {'por022a08_p_SPN'}, {'por022a12_p_SPN'});

global STACK XXX MODULES NARF_PATH NARF_SAVED_MODELS_PATH ...
    NARF_MODULES_PATH NARF_SAVED_ANALYSIS_PATH;


MODULES = scan_directory_for_modules(NARF_MODULES_PATH);

mm = {}; 

% GROUP 1: MODEL LOADING
mm{1} = [];
mm{1}.env200hz = {MODULES.load_stim_resps_from_baphy.mdl(...
                                    struct('raw_resp_fs', 200, ...
                                           'raw_stim_fs', 200,...
                                           'stimulus_format','envelope'))};
mm{1}.env100hz = {MODULES.load_stim_resps_from_baphy.mdl(...
                                    struct('raw_resp_fs', 100, ...
                                           'raw_stim_fs', 100,...
                                           'stimulus_format','envelope'))};                                    
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
%                                        
% GROUP 2: BEFORE THE LINEAR FILTER
mm{2} = [];
mm{2}.none = {MODULES.passthru};
mm{2}.sqrt = {MODULES.nonlinearity.mdl(struct('phi', [], ...
                                              'nlfn', @(phi, z) sqrt(z)))};
mm{2}.log  = {MODULES.nonlinearity.mdl(struct('phi', [], ...
                                              'nlfn', @(phi, z) log(z + 10^-6)))};
mm{2}.volterra = {MODULES.concat_second_order_terms}; 
mm{2}.depress  = {MODULES.depression_filter_bank}; 

% GROUP 3: THE FIR FILTER 
mm{3} = [];
mm{3}.fir =      {MODULES.fir_filter.mdl(struct('num_coefs', 12))};

mm{3}.norm_fir = {MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12))};
              
mm{3}.inhibexcit = {MODULES.normalize_channels, ...
                    MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                  'output', 'inhib')), ...
                    MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                  'output', 'excit')), ...                                                  
                    MODULES.nonlinearity.mdl(struct('input_stim', 'inhib', ...
                                                    'output', 'inhib', ...
                                                    'phi', [0], ...
                                                    'nlfn', @(phi, z) - zero_below_thresh(z))), ...
                    MODULES.nonlinearity.mdl(struct('input_stim', 'excit', ...
                                                    'output', 'excit', ...
                                                    'phi', [0], ...
                                                    'nlfn', @zero_below_thresh)), ...
                    MODULES.sum_fields.mdl(struct('inputs', {{'inhib', 'excit'}}, ...
                                                  'output', 'stim'))};

% GROUP 4: POST-FILTER NONLINEARITIES
mm{4} = [];
mm{4}.none = {MODULES.passthru};
mm{4}.exp  = {MODULES.nonlinearity.mdl(struct('phi', [1 1], ...
                                              'nlfn', @exponential))};
mm{4}.sig  = {MODULES.nonlinearity.mdl(struct('phi', [0 1 1 1], ...
                                              'nlfn', @sigmoidal))};
mm{4}.step = {MODULES.nonlinearity.mdl(struct('phi', [0], ...
                                              'nlfn', @zero_below_thresh))};
mm{4}.poly = {MODULES.nonlinearity.mdl(struct('phi', [1 1 1 1], ...
                                              'nlfn', @polyval))};       

% GROUP 5: PERFORMANCE METRICS TO OPTIMIZE ON
% Note that you should include all performance metrics here but choose only
% one to be used as the objective function
mm{5} = [];
mm{5}.mse  = {MODULES.correlation, ...
              MODULES.mean_squared_error.mdl(struct('output', 'score'))};
mm{5}.corr = {MODULES.correlation.mdl(struct('output', 'score')), ...
              MODULES.mean_squared_error};

% ------------------------------------------------------------------------
% BUILD THE MODELS

opts = cellfun(@fieldnames, mm, 'UniformOutput', false);
N_opts = cellfun(@(m) length(fieldnames(m)), mm);
N_models = prod(N_opts);
fprintf('Number of models to be tested: %d\n', N_models); 

results = {};

for ii = 1:N_models,
    STACK = {};
    XXX = {};
    XXX{1}.cellid = cellid;
    XXX{1}.training_set = training_set;
    XXX{1}.test_set = test_set;

    exit_codes = [];
    tic;    
    
    [i1 i2 i3 i4 i5 i6 i7 i8 i9] = ind2sub(N_opts, ii);
    indexes = [i1 i2 i3 i4 i5 i6 i7 i8 i9];
    opt_names = arrayfun(@(gi, ind) opts{gi}{ind}, ...
                      1:length(opts), indexes(1:length(opts)), ...
                     'UniformOutput', false);

    blocks = cellfun(@getfield, mm, opt_names, 'UniformOutput', false);
   
    tmp = cellfun(@(n) sprintf('%s_', n), opt_names, 'UniformOutput', false);
	modelname = strcat(tmp{:}, cellid);
    fprintf('\n\nMODEL [%d/%d]: %s\n', ii, N_models, modelname);
    
     for jj = 1:length(blocks),
         fprintf('Fitting stage [%d/%d]\n', jj, length(blocks));
         STACK = cat(1, STACK, blocks{jj}');
         if jj < 3
             exit_codes(jj) = NaN;
         else
             recalc_xxx(1);
%             exit_codes(jj) = fit_with_lsqcurvefit();
         end
     end
%     
%     results{ii}.score_train_corr = XXX{end}.score_train_corr;
%     results{ii}.score_test_corr = XXX{end}.score_test_corr;
%     results{ii}.score_train_mse = XXX{end}.score_train_mse;
%     results{ii}.score_test_mse = XXX{end}.score_test_mse;
%     firmod = find_module(STACK, 'fir_filter');
%     results{ii}.fir_coefs = firmod.coefs;
%     
%     % TODO: Find a 'nonlinearity' module in the stack, but start the search
%     % at a lower depth than just 1 so that it doesn't find the pre-FIR
%     % nonlinearity. 
%     % results{ii}.nl_phi = STACK{find_module(STACK{3}'nonlinearity')}.phi;
%     results{ii}.fit_time = toc;
%     results{ii}.n_free_params = length(pack_fittables(STACK));
%     results{ii}.exit_codes = exit_codes;
%     
%     results{ii}.cellid = XXX{1}.cellid;
%     results{ii}.training_set = XXX{1}.training_set;
%     results{ii}.test_set = XXX{1}.test_set;
%     results{ii}.optimized_with = 'lsqcurvefit';
%     
%     save_model_stack([NARF_SAVED_MODELS_PATH filesep modelname '.mat'], STACK, XXX);
end

% Save the results matrix to file
save([NARF_SAVED_ANALYSIS_PATH filesep cellid '_results.mat'], 'results');
