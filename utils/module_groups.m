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
            mm = {MODULES.nonlinearity.mdl(struct('phi', [2], 'nlfn', @rootcompressor))};
        case {'root3'} 
            mm = {MODULES.nonlinearity.mdl(struct('phi', [3], 'nlfn', @rootcompressor))};
        case {'root4'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [4], 'nlfn', @rootcompressor))};
        case {'root5'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [5], 'nlfn', @rootcompressor))};
        case {'rootfit'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [2], 'nlfn', @rootcompressor))};
                                             
        case {'log1'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-1], 'nlfn', @logcompressor))};
        case {'log2'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-2], 'nlfn', @logcompressor))};
        case {'log3'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-3], 'nlfn', @logcompressor))};
        case {'log4'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-4], 'nlfn', @logcompressor))};
        case {'log5'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-5], 'nlfn', @logcompressor))};
        case {'logfit'}
            mm = {MODULES.nonlinearity.mdl(struct('phi', [-5], 'nlfn', @logcompressor, ...
                                                  'fit_fields', {{'phi'}}))};

        case {'fir'} 
            mm = {MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs'}}))};
        case {'firb'} 
            mm = {MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                                'fit_fields', {{'coefs', 'baseline'}}))};
        case {'depfir'}
            mm = {MODULES.depression_filter_bank, ...
                  MODULES.normalize_channels, ...
                  MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                               'fit_fields', {{'coefs', 'baseline'}}))};                                          
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
                                                   'nlfn', @(phi, z) - zero_below_thresh(phi, z))), ...
                   MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                   'input_stim', 'excit', ...
                                                   'output', 'excit', ...
                                                   'phi', [0], ...
                                                   'nlfn', @zero_below_thresh)), ...
                   MODULES.sum_fields.mdl(struct('inputs', {{'inhib', 'excit'}}, ...
                                                 'output', 'stim'))};
                                             
        case {'nonl'}
            mm = {MODULES.passthru};
        case {'npnl'}
            mm = {MODULES.nonparm_nonlinearity};
        case {'sig'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [0 1 1 0], ...
                                                  'nlfn', @sigmoidal))};
        case {'exp'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [1 1], ...
                                                  'nlfn', @exponential))};
        case {'step'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [0], ...
                                              'nlfn', @zero_below_thresh))};
        case {'poly'}
            mm = {MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [1 1 1 1], ...
                                              'nlfn', @polyval))};       

        case {'mse'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score'))};                     
        case {'msereg'}
            mm = {MODULES.mean_squared_error.mdl(struct('output', 'score', ...
                                                        'smoothness_weight', 10^-6))};
            
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
        otherwise
            error('WTF kind of key is that: %s', g);
    end
    
    ret.(g) = mm;
    
end
             
