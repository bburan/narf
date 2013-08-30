function cp1()

global MODULES STACK META;

append_module(MODULES.nonlinearity.mdl(struct('phi', [1 10 -2], ...
                                              'nlfn', @nl_linexp, ...
                                              'fit_fields', {{'phi'}})));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 1, 'fit_fields', {{}})));
append_module(MODULES.mean_squared_error);
META.perf_metric = @pm_mse;

% There may be multiple channels, so do the correlation with each
ll = size(STACK{end-1}{1}.coefs, 1);
thebst = nan; 
thephi = pack_fittables(STACK);
for ii = 1:ll
    STACK{end-1}{1}.coefs = zeros(ll, 1);
    STACK{end-1}{1}.coefs(ii, 1) = 1;
    calc_xxx(2);
    % fit_fminlsq();
    fit_fminsearch(optimset('MaxIter', 1000, 'MaxFunEvals', 1000, ...
                            'TolFun', 1e-12, 'TolX', 1e-9));    
    if isnan(thebst) || thebst < META.perf_metric()
        thebst = META.perf_metric();
        thephi = pack_fittables(STACK);
    end
end
unpack_fittables(thephi);
pop_module();
pop_module();
STACK{end}{1} = rmfield(STACK{end}{1}, 'fit_fields');