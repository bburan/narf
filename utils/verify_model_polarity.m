function verify_model_polarity()
% No arguments. 
% Flips the sign of every element in the FIR filter iff the corresponding
% NPNL, NPNLX, NPFNL, SENL, or GMM nonlinearity output has a negative slope
% and the NL module is two modules after the FIR filter

global STACK XXX META;

[firmod, firmod_idx] = find_modules(STACK, 'fir_filter', true);

firmod.plot_fns{1}.fn(STACK(1:firmod_idx), XXX(1:firmod_idx+1));

nlidx = firmod_idx + 2; % TODO: This assumes firn or depn, which is not always true!

if strcmp(STACK{nlidx}.name, 'nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'sparse_empirical_nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'gmm_nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'nonparm_filter_nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'nonparm_nonlinearity_x') | ...
        strcmp(STACK{nlidx}.name, 'nonparm_nonlinearity')
        
    x = XXX{nlidx};
    mdl = STACK{nlidx};
    
    % Flip the sign of the coefs if there is a negative slope
    V1 = [];
    V2 = [];
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};
        V1 = cat(1, V1, x.dat.(sf).(mdl.input_stim)(:));
        V2 = cat(1, V2, x.dat.(sf).(mdl.input_resp)(:));
    end
    D = [V1(:) V2(:)]; 
    D = sortrows(D);
    P = polyfit(D(:,1),D(:,2), 1);
    % R = corrcoef(excise([V1 V2]));
    if (P(1) < 0)
        STACK{firmod_idx}.coefs = - STACK{firmod_idx}.coefs;
    end
    recalc_xxx(firmod_idx); 
    
end
