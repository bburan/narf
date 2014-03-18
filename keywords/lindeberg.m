function lindeberg()
% March 2014 - lienard
% modified from 'fir.m'

% for now, this is mostly a test to try to get to a symmetrical STRF
% -- not a real implementation of Linderberg's family of functions

global MODULES XXX;
global ANESTHETIZED_MOUSE

if ANESTHETIZED_MOUSE==1,
    disp('fir: Using longer num_coefs for mouse data.');
    dep_tau_norm=100;
    fir_num_coefs=20;
    stop_exp=3.0;
else
    dep_tau_norm=100;
    fir_num_coefs=20; % was 12 previously (JL)
    stop_exp=2.0;
end

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.lindeberg_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'baseline',meanresp,...    
                                'fit_fields', {{'lincoefs','baseline'}})));
%                                 'fit_fields', {{'lincoefs'}})));

fitSubstack([],10^-(stop_exp));
