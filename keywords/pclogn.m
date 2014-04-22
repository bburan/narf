function pclogn()

global MODULES STACK XXX;

mdl = STACK{end}{1}
x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, C] = size(x.dat.(sf).(mdl.output)); 

append_module(MODULES.nonlinearity.mdl(struct('phi', repmat([-2], C, 1), ...
                                              'nlfn', @nl_log, ...
                                              'fit_fields', {{'phi'}})));
STACK{end}{1}.auto_plot=[];
