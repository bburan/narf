function firtemp()
% Fit a FIR filter temporarily
global STACK MODULES XXX;
signal = 'stim';
fns = fieldnames(XXX{end}.dat);

fir_num_coefs=10;
n_input_chans = size(XXX{end}.dat.(fns{1}).(signal), 3);
append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'baseline',meanresp,...
                                'fit_fields', {{'coefs','baseline'}})));
fitSubstack(length(STACK), 10^-4);