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

% if there are zero or multiple filters, this is a more complex-form model
% and not to be messed with here.
if length(firmod_idxs) > 1 || isempty(firmod_idxs)
    return
end

% If there are any free parameters after the FIR filter, we should NOT
% try to flip anything, because it would mess up the model.
for kk = (1+firmod_idxs{1}):length(STACK)
    for ii = 1:length(STACK{kk})
        if isfield(STACK{kk}{ii}, 'fit_fields') && ~isempty(STACK{kk}{ii}.fit_fields)
            fprintf('Free parameters found after the FIR filter. Flipping the polarity will probably break the model, so I am refusing!\n');
            return;
        end
    end
end

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
   
    [~, idxs, ~] = unique(V1, 'first'); % Remove duplicates
    D = [V1(idxs) V2(idxs)]; 
    D = excise(D);
    D = sortrows(D);
    winsize = ceil(size(D, 1)/100);
    if winsize < 1
        fprintf('RESPAVG and fir filter outputs are all NaN, skipping\n');
        return;
    end
    D = conv_fn(D, 1, @nanmean, winsize, 0);
    xs = D(:,1);
    ys = D(:,2);
    
    if length(ys) == 1
        fprintf('verify_model_polarity.m: All values equal. Skipping.\n');
        return;
    end
    
    idx10 = ceil(length(ys) / 20);
    y1 = mean(ys(1:idx10));
    y2 = mean(ys(end-idx10:end));
    if (y1 > y2)
        fprintf('verify_model_polarity.m: Detected negative slope. Flipping.\n');
  
        % There are no later fittable params, we can flip safely
        for aa = 1:length(STACK{idx})
            STACK{idx}{aa}.coefs = - STACK{idx}{aa}.coefs;
        end
    end
        
    calc_xxx(idx); 
    
end
