function initrnd()
% Initializes all FIR filters modules with the 'coefs' tag in 'fit_fields'
% such that the coefs for that module are all between 10^-3 and zero. 
% This is essentially the same as init0, but with nonzero starting values
% to help avoid all-zero fitting bugs.

global STACK;

[~, mod_idxs] = find_modules(STACK, 'fir_filter');

for ii = 1:length(mod_idxs)
    for jj = 1:length(STACK{mod_idxs{ii}})
        if isfield(STACK{mod_idxs{ii}}{jj}, 'fit_fields') && ...
           any(strcmp('coefs', STACK{mod_idxs{ii}}{jj}.fit_fields))
            STACK{mod_idxs{ii}}{jj}.coefs = normrnd(0, 10^-3, size(STACK{mod_idxs{ii}}{jj}.coefs));
        end
    end
end
