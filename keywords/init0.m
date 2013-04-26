function init0()
% Initializes all FIR filters modules with the 'coefs' tag in 'fit_fields'
% such that the coefs for that module are all zeros.

global STACK;

[~, mod_idxs] = find_modules(STACK, 'fir_filter');

% Initialize coefs to all ones 
for ii = 1:length(mod_idxs)
    for jj = 1:length(STACK{mod_idxs{ii}})
        if isfield(STACK{mod_idxs{ii}}{jj}, 'fit_fields') && ...
           any(strcmp('coefs', STACK{mod_idxs{ii}}{jj}.fit_fields))
            STACK{mod_idxs{ii}}{jj}.coefs = zeros(size(STACK{mod_idxs{ii}}{jj}.coefs));
        end
    end
end
