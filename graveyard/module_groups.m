function ret = module_groups(varargin)
% Returns a struct with modules under each field that correspond to the
% uniquely identified "group name". Read on if that doesn't make sense.
%
% MOTIVATION:
%     Ninety percent of the time, we will want the same sorts of module
%     groupings to be used in different analyses. For example, most of the
%     time we want a depression bank filter to be followed by a normalizing
%     filter and then a FIR filter. If we call this 'depfir' as a
%     shorthand, we can refer to this grouping of modules in a shorthand
%     way. This is convenient in many cases. 
%     
%     Shorthand notation is good, but if we end up cutting and pasting code
%     from one analysis file to another, our code can still turn into a big
%     dangerous ball. If we changed the depression time constants and r
%     removed the normalization, yet still refer to it simply as 'depfir', 
%     we would be mislead when comparing the performance between another
%     analysis for which the name 'depfir' means something else.
%
%     Also, NARF uses these shorthand group names during file naming,
%     analysis caching, and other conveniences. 
% 
%     Therefore, this file contains all the definitions of common module 
%     groups in a single place, so that a unique identifier stays unique. 
%     It is intended that as certain groupings of modules become used by 
%     more than one analysis file, you will edit this file and add a new
%     identifier to this list, leaving the old ones in place. 
%    
%     Generally speaking, if you CHANGE the definition of any of these, you
%     should DELETE ALL CORRESPONDING MODEL FILES so that they can be
%     automatically rebuilt and your shorthand names stay in sync with the
%     saved model files.
% 
% USE:
%     GROUP_NAMES    A cell array of strings (which identify module
%                    groupings uniquely).   
%
% TODO: At a later date, perhaps it would be better to break off all these
% definitions into a user-controlled file or directory so that new module
% groups could be defined without needing to edit this file. For now, the
% convenience of keeping group names here wins out, however.
%
% TODO: Also, right now underscores are NOT ALLOWED in shorthand names.

global MODULES;

if isempty(MODULES)
    MODULES = scan_directory_for_modules(NARF_MODULES_PATH);
end

group_names = varargin;

ret = {};

for ii = 1:length(group_names)
    g = group_names{ii};
    
    switch g
        case {'env100'}
            mm = {MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 100, ...
                                  'raw_stim_fs', 100,...
                                  'stimulus_format', 'envelope'))}; 
        case {'env100n'}
            mm = {MODULES.load_stim_resps_from_baphy.mdl(...
                           struct('raw_resp_fs', 100, ...
                                  'raw_stim_fs', 100,...
                                  'include_prestim', false,...
                                  'stimulus_format', 'envelope'))}; 
        case {'env200'}
            mm = {MODULES.load_stim_resps_from_baphy.mdl(...
                                     struct('raw_resp_fs', 200, ...
                                            'raw_stim_fs', 200,...
                                            'stimulus_format', 'envelope'))};     
        
        case {'ellip100'}
            mm =  {MODULES.load_stim_resps_from_baphy.mdl(...
                                    struct('raw_resp_fs', 200, ...
                                           'raw_stim_fs', 100000)), ...
                  MODULES.elliptic_bandpass_filter_bank, ...
                  MODULES.downsample_with_fn.mdl(...
                                    struct('downsampled_freq', 200, ...
                                           'conv_fn', @mean, ...
                                           'postconv_fn', @abs))};
        case {'ellip200'}
            mm = {MODULES.load_stim_resps_from_baphy.mdl(...
                                     struct('raw_resp_fs', 200, ...
                                            'raw_stim_fs', 100000)), ...
                   MODULES.elliptic_bandpass_filter_bank, ...
                   MODULES.downsample_with_fn.mdl(...
                                     struct('downsampled_freq', 200, ...
                                            'postconv_fn', @abs))};

        case {'nocomp'}
            mm = {MODULES.passthru};
    
        case {'root2'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [2], 'nlfn', @nl_root))};
        case {'root3'} 
            mm = {MODULES.nonlinearity.mdl(struct('phi', [3], 'nlfn', @nl_root))};
        case {'root4'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [4], 'nlfn', @nl_root))};
        case {'root5'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [5], 'nlfn', @nl_root))};
        case {'rootfit'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [2], 'nlfn', @nl_root))};
                                             
        case {'log1'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-1], 'nlfn', @nl_log))};
        case {'log2'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-2], 'nlfn', @nl_log))};
        case {'log3'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-3], 'nlfn', @nl_log))};
        case {'log4'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-4], 'nlfn', @nl_log))};
        case {'log5'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-5], 'nlfn', @nl_log))};
        case {'logfit'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-5], 'nlfn', @nl_log, ...
                                                  'fit_fields', {{'phi'}}))};
        case {'log1b'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-1 -log(0+10^-1)], 'nlfn', @nl_log))};
        case {'log2b'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-2 -log(0+10^-2)], 'nlfn', @nl_log))};
        case {'log3b'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-3 -log(0+10^-3)], 'nlfn', @nl_log))};          
        case {'log4b'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-4 -log(0+10^-4)], 'nlfn', @nl_log))};
        case {'log5b'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-5 -log(0+10^-5)], 'nlfn', @nl_log))};

        case {'l1'}
            mm = {MODULES.smooth_respavg, ...
                  MODULES.nonlinearity.mdl(struct('phi', [-1 -log(0+10^-1)], 'nlfn', @nl_log))};
        case {'l2'}
            mm = {MODULES.smooth_respavg, ...
                  MODULES.nonlinearity.mdl(struct('phi', [-2 -log(0+10^-2)], 'nlfn', @nl_log))};
        case {'l3'}
            mm = {MODULES.smooth_respavg, ...
                  MODULES.nonlinearity.mdl(struct('phi', [-3 -log(0+10^-3)], 'nlfn', @nl_log))};          
        case {'l4'}
            mm = {MODULES.smooth_respavg, ...
                  MODULES.nonlinearity.mdl(struct('phi', [-4 -log(0+10^-4)], 'nlfn', @nl_log))};
        case {'l5'}
            mm = {MODULES.smooth_respavg, ...
                  MODULES.nonlinearity.mdl(struct('phi', [-5 -log(0+10^-5)], 'nlfn', @nl_log))};

        case {'fir'} 
            mm = {MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs'}}))};
        case {'firb'} 
            mm = {MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs', 'baseline'}}))};
        case {'firc'} 
            mm = {MODULES.normalize_channels.mdl(struct('force_positive', true)), ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs', 'baseline'}}))};
        case {'firn'}  % Normalized, no baseline
            mm = {MODULES.normalize_channels.mdl(struct('force_positive', true)), ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs'}})), ...
                  MODULES.normalize_channels};
            
                                            
        case {'depfir'}
            mm = {MODULES.depression_filter_bank, ...
                  MODULES.normalize_channels.mdl(struct('force_positive', true)), ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                               'fit_fields', {{'coefs', 'baseline'}}))};                                          
                                                  
      case {'depn'}
            mm = {MODULES.normalize_channels.mdl(struct('force_positive', true)), ...
                  MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0,  0.2   0.2    0.2    0.2], ...
                           'tau',      [0,  25,   75,    150,   400], ...
                           'pedestal',0.0,...
                           'num_channels', 5)), ...
                  MODULES.normalize_channels.mdl(struct('force_positive', true)), ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 10, ...
                                               'fit_fields', {{'coefs'}})), ...
                  MODULES.normalize_channels.mdl(struct('force_positive',true))};
      
      case {'depped'}
            mm = {MODULES.normalize_channels.mdl(struct('force_positive', true)), ...
                  MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0,  0.2   0.2    0.2    0.2], ...
                           'tau',      [0,  25,   75,    150,   400], ...
                           'pedestal',0.1,...
                           'num_channels', 5)), ...
                  MODULES.normalize_channels.mdl(struct('force_positive', true)), ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 10, ...
                                               'fit_fields', {{'coefs'}})), ...
                  MODULES.normalize_channels.mdl(struct('force_positive',true))};
            
        case {'voltfir'} 
            mm = {MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs'}}))};                                            
                                            
        case {'inex'}
            mm  = {MODULES.normalize_channels, ...
                   MODULES.fir_filter.mdl(struct('fit_fields', {{'coefs'}}, ...
                                                 'num_coefs', 12, ...
                                                 'output', 'inhib')), ...
                   MODULES.fir_filter.mdl(struct('fit_fields', {{'coefs'}}, ...
                                                 'num_coefs', 12, ...
                                                 'output', 'excit')), ...                                                  
                   MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                   'input_stim', 'inhib', ...
                                                   'output', 'inhib', ...
                                                   'phi', [0], ...
                                                   'nlfn', @(phi, z) - nl_zerothresh(phi, z))), ...
                   MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                   'input_stim', 'excit', ...
                                                   'output', 'excit', ...
                                                   'phi', [0], ...
                                                   'nlfn', @nl_zerothresh)), ...
                   MODULES.sum_fields.mdl(struct('inputs', {{'inhib', 'excit'}}, ...
                                                 'output', 'stim'))};
                                             
        case {'nonl'}
            mm = {MODULES.passthru};
        case {'npnl'}
            mm = {MODULES.nonparm_nonlinearity};
        case {'npnlx'}
            mm = {MODULES.nonparm_nonlinearity_x};
        case {'npfnl'}
            mm = {MODULES.nonparm_filter_nonlinearity};
        case {'npfnl3'}
            mm = {MODULES.nonparm_filter_nonlinearity.mdl(struct('gwinval', 3))};
        case {'npfnl4'}
            mm = {MODULES.nonparm_filter_nonlinearity.mdl(struct('gwinval', 4))};
        case {'npfnl5'}
            mm = {MODULES.nonparm_filter_nonlinearity.mdl(struct('gwinval', 5))};
        case {'senl'}
            mm = {MODULES.sparse_empirical_nonlinearity};
        case {'se2d'}
            mm = {MODULES.sparse_empirical_nonlinearity_2d};
        case {'senl25'}
            mm = {MODULES.sparse_empirical_nonlinearity.mdl(struct('relvar', 0.25))};
        case {'senl3'}
            mm = {MODULES.sparse_empirical_nonlinearity.mdl(struct('relvar', 0.3))};

        case {'gmm3'}
            mm = {MODULES.gmm_nonlinearity.mdl(struct('num_pts', 500, ...
                                                      'num_gaussians', 3))};
        case {'gmm4'}
            mm = {MODULES.gmm_nonlinearity.mdl(struct('num_pts', 500, ...
                                                      'num_gaussians', 4))};            
        case {'gmm5'}
            mm = {MODULES.gmm_nonlinearity.mdl(struct('num_pts', 500, ...
                                                      'num_gaussians', 5))};
        case {'gmm6'}
            mm = {MODULES.gmm_nonlinearity.mdl(struct('num_pts', 500, ...
                                                      'num_gaussians', 6))};
            
        case {'sig'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [0 1 1 0], ...
                                                  'nlfn', @nl_sigmoid))};
        case {'exp'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [1 1], ...
                                                  'nlfn', @nl_exponential))};
        case {'step'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [0], ...
                                              'nlfn', @nl_zerothresh))};
        case {'poly'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [1 1 1 1], ...
                                              'nlfn', @polyval))};       

        case {'mse'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score'))};                     
        case {'msereg'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'smoothness_weight', 10^-6))};
        case {'mses0'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'sparseness_weight', 10^0))};
        case {'mses1'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'sparseness_weight', 10^-1))};
                                                    
        case {'mses2'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'sparseness_weight', 10^-2))};
                                                    
        case {'mses3'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'sparseness_weight', 10^-3))};     
        
        case {'mses4'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'sparseness_weight', 10^-4))};
        case {'mses5'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'sparseness_weight', 10^-5))};
        case {'mses6'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'sparseness_weight', 10^-6))};
        case {'fmin'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_fminsearch))};    
        case {'lsq'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_lsq))};
        case {'fminlsq'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_fminlsq))};   
        case {'boost'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_boost))};
        case {'lsqnl'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_lsqnonlin))};
        case {'twostep'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_twostep))};
        case {'fminunc'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_fminunc))};
        case {'slsq'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_slsq))};
        case {'sboost'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_sboost))};
        case {'slsqtwo'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_slsqtwo))};
        case {'anneal'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_anneal))};
        case {'gene'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_genetic))};
        case {'sgene'}
            mm = {MODULES.correlation.mdl(struct('fitter', @fit_sgene))};
            
        case {'nlsep1'}
            mm = {MODULES.correlation.mdl(struct('fitter', @( )fit_nlsep(1)))};
        case {'nlsep2'}
            mm = {MODULES.correlation.mdl(struct('fitter', @() fit_nlsep(2)))};
        case {'nlsep3'}
            mm = {MODULES.correlation.mdl(struct('fitter', @() fit_nlsep(3)))};
        case {'nlsep4'}
            mm = {MODULES.correlation.mdl(struct('fitter', @() fit_nlsep(4)))};
        case {'nlsep5'}
            mm = {MODULES.correlation.mdl(struct('fitter', @() fit_nlsep(5)))};
        case {'sb'}
            mm = {MODULES.correlation.mdl(struct('fitter', @() fit_sparsebayes))};
        otherwise
            error('WTF kind of key is that: %s', g);
    end
    
    ret.(g) = mm;
    
end
             
