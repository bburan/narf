function verify_model_polarity()
% verify_model_polarity()
%
% (Requires that the STACK and XXX data structures be fully initialized.)
%
% For each FIR filter in the STACK, 
%     Flips the sign of every element in that filter iff its output is
%     creating a scatter plot with a negative slope.
%
% No arguments or return values.

global STACK XXX;

[~, firmod_idxs] = find_modules(STACK, 'fir_filter');

for kk = 1:length(firmod_idxs)
    idx = firmod_idxs{kk};
    
    firmod = STACK{idx}{1};
    x = XXX{idx+1};
        
    % Flip the sign of the coefs if there is a negative slope
    V1 = [];
    V2 = [];
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};
        V1 = cat(1, V1, x.dat.(sf).(firmod.output)(:));
        V2 = cat(1, V2, x.dat.(sf).respavg(:));
    end
    D = [V1(:) V2(:)]; 
    D = sortrows(D);
    P = polyfit(D(:,1),D(:,2), 1);
    
    if (P(1) < 0)
        for aa = 1:length(STACK{idx})
            STACK{idx}{aa}.coefs = - STACK{idx}{aa}.coefs;
        end
    end
    calc_xxx(idx); 
    
end
