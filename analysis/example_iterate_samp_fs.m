% A script to demonstrate how to use NARF and 'narf_modelpane.m'
%
% Essentially, this script sets up a model using a gammatone filter bank
% prefilter, a simple decimation using the mean value and square root, and
% an FIR filter which accepts the output of every gammatone filter. When
% you feed it some data, boosting is used to estimate the parameters of
% a FIR filter. The FIR filter coefficients should then define an STRF.
%
% Note that this is a very primitive first attempt, and should not be used
% for any real scientific use, mostly because the boosting routine is
% really really simple and naive.

narf_set_path;

global NARF_PATH NARF_SAVED_MODELS_PATH STACK XXX;

% Build a list of the files in the 'modules/' directory
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules/']);

% XXX defines the initial data
XXX = {};
XXX{1} = [];
XXX{1}.cellid = 'por022b-a1';
XXX{1}.training_set = {'por022b02_p_TOR'};
XXX{1}.test_set = {};

% For several different sampling rates, 
sampling_rates = {100, 150, 200, 250};

for sr = sampling_rates, sr = sr{1};
    
    STACK{1} = mdls.load_stim_resps_from_baphy.mdl(struct('raw_resp_fs', sr));
    STACK{2} = mdls.gammatone_filter_bank.mdl(struct('num_channels', 10, ...
                                                     'bank_min_freq', 250, ...
                                                     'bank_max_freq', 8000));
    STACK{3} = mdls.downsample_with_fn.mdl(struct('downsampled_freq', sr));
    STACK{4} = mdls.fir_filter.mdl(struct('num_dims', 10, ...
                                          'fit_fields', {{'coefs'}}));
                                      
    recalc_xxx(1);

    %%%%%%%%%%%%%%%%%%%%%%%%5
    % UNCOMMENT TO USE IVAR'S FIT
    %     
    %     n_iters = 50;
    %     x_0 = pack_fittables(STACK);
    %     if isempty(x_0)
    %         log_msg('No parameters were selected to be fit.');
    %         return;
    %     end
    %     
    %     [x_bst, s_bst] = boosting(x_0', @correlation_of_downsampled_signals, ...
    %                              @(n,x,s)(n > n_iters), 1.0);
    %     
    %     unpack_fittables(x_bst);
    %     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UNCOMMENT TO USE STEPHENS FIT
    
    m = [];
    m.num_dims = 2;
    m.altcore = 'cdcore';
    m.maxlag = 19;
    m.resampcount = 19;  % Why is the total number of things this + 1?
    m.sfscount = 10;
    m.sfsstep = 3;
    m.rasterfs = sr; 
    m.coefs = zeros(m.num_dims, 10);
    
    % Run Stephen's fitting routine at depth 4
    STACK{4}.coefs = do_stephen_fit(m, XXX{4});
    ss = size(STACK{4}.coefs);
    STACK{4}.num_dims = ss(1);
    STACK{4}.num_coefs = ss(2);
    % Recompute now just the last little bit
    recalc_xxx(4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Save the model struct
    filename = [NARF_SAVED_MODELS_PATH filesep 'svd_sr_' num2str(sr) '.mat']
    
    save_model_stack(filename, STACK);
    
end
