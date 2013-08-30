function inita()
% Initializes all FIR filters modules with the 'coefs' tag in 'fit_fields'
% such that the coefs for that module are all zeros, except for
% coefficients 2,3,4 of each channel, which are a 1.

global STACK;

[~, mod_idxs] = find_modules(STACK, 'fir_filter');

for ii = 1:length(mod_idxs)
    for jj = 1:length(STACK{mod_idxs{ii}})
        if isfield(STACK{mod_idxs{ii}}{jj}, 'fit_fields') && ...
           any(strcmp('coefs', STACK{mod_idxs{ii}}{jj}.fit_fields))
            STACK{mod_idxs{ii}}{jj}.coefs = zeros(size(STACK{mod_idxs{ii}}{jj}.coefs));
            STACK{mod_idxs{ii}}{jj}.coefs(:,2) = 1;
            STACK{mod_idxs{ii}}{jj}.coefs(:,3) = 1;
            STACK{mod_idxs{ii}}{jj}.coefs(:,4) = 1;
        end
    end
end
