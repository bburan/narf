function lindeberg()
% March 2014 - lienard
% modified from 'fir.m'

% NOT TO BE USED %
% kept as a way to test different hypotheses locally

global MODULES XXX;

fir_num_coefs=20;
stop_exp=2.0;


append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.lindeberg_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'order_x',0,...
                                'order_t',0,...
                                'baseline',meanresp,...    
                                'fit_fields', {{'lincoefs'}})));
%                                 'fit_fields', {{'lincoefs','baseline'}})));

fitSubstack([],10^-(stop_exp));
