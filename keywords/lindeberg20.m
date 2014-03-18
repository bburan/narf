function lindeberg20()
% March 2014 - lienard
% modified from 'fir.m'

% Linderberg's time-causal kernels of order (dx=2s,dt=0)

global MODULES XXX;

dep_tau_norm=100;
fir_num_coefs=20; % was 12 previously (JL)
stop_exp=2.0;


append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.lindeberg_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'order_x',2,...
                                'order_t',0,...
                                'baseline',meanresp,...
                                'fit_fields', {{'lincoefs','baseline'}})));
%                                 'fit_fields', {{'lincoefs'}})));

fitSubstack([],10^-(stop_exp));
