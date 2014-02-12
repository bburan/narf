% zthresh3- same as zthresh2 but no fitSubstack 
function zthresh3()

global MODULES XXX;

meanstim = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'stim'));
stdstim = nanstd(flatten_field(XXX{end}.dat,XXX{end}.training_set,'stim'));
mm=meanstim-stdstim*2;
meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
stdresp = nanstd(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
rr=meanresp-stdresp*2;

append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [mm 100 rr], ...
                                              'nlfn', @nl_zerothresh)));
%fitSubstack();