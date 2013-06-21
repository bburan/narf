function logfit()

global MODULES STACK;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 -log(0+10^-2)], ...
                                              'nlfn', @nl_log, ...
                                              'fit_fields', {{'phi'}})));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 1)));
append_module(MODULES.correlation);
META.perf_metric = @pm_corr;
% There may be multiple channels, so do the correlation with each
ll = size(STACK{end-2}{1}.coefs, 2);
themax = 0; 
thephi = pack_fittables();
for ii = 1:ll
    STACK{end-2}{1}.coefs = zeros(1, ll);
    STACK{end-2}{1}.coefs(1, ii) = 1;
    fit_fminlsq();
    if themax < META.perf_metric()
        themax = META.perf_metric();
        thephi = pack_fittables();
    end
end
unpack_fittables(thephi);
pop_module();
pop_module();
STACK{end}{1} = rmfield(STACK{end}{1}, 'fit_fields');