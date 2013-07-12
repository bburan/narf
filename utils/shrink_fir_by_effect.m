function shrink_fir_by_effect(N_fraction)
% Ad hoc hack that trades N_FRACTION of performance metric quality for
% sparsity. 
%
% Shrinks coefficients of all FIR filter in all parameter sets towards zero
% based upon how much a small perturbation in them affects the performance
% metric. 

if ~exist('N_fraction', 'var')
    N_fraction = 0.001;
end

global STACK META;

[~, mod_idxs] = find_modules(STACK, 'fir_filter');

coefs={};
effects = {};
init_perf = META.perf_metric(); 
rms_pwr = {};
stack_cache = STACK;
for ii = 1:length(mod_idxs)
    idx = mod_idxs{ii};
    for jj = 1:length(STACK{idx})
        coefs{ii}{jj} = STACK{idx}{jj}.coefs;
        effects{ii}{jj} = ones(size(STACK{idx}{jj}.coefs));
        
        for kk = 1:length(coefs{ii}{jj}(:))
            % Perturb coefficient by 1%
            STACK{idx}{jj}.coefs(kk) = 1.01 * STACK{idx}{jj}.coefs(kk);
            calc_xxx(idx);
            
            % Estimate the relative effect on the output and save it
            effects{ii}{jj}(kk) = abs(init_perf - META.perf_metric()) / init_perf; 
            
            % Restore the old coefficient
            STACK{idx}{jj}.coefs(kk) = coefs{ii}{jj}(kk);
        end
        
        effects{ii}{jj}(:) = effects{ii}{jj}(:) ./ mean(effects{ii}{jj}(:));
    end
end

% Three options to modify the coefficients based on their effects:
% 1. Discrete chopping: Zero out bottom n% most ineffectual coefficients
% 2. Scaled: Multiply everything by its effectivity, then rescale
% 3. Nonlinear: Rescale According to its rank, sqrt of the effect, etc
%
% Difficult Problem with Option #2:
% If I scale all the coefficients and don't have a normalization module
% between the FIR filter and a parametric nonlinearity, the previously fit
% nonlinearity will no longer have the appropriate domain/range to
% transform the FIR output into something useful. 
% 
% I'm going to therefore punt on this, and assume that we are just doing
% this for FIRN (normalized FIR on both sides) instead of FIRNO (just
% normalization on the input side).
%
% In some respects, just zeroing out the unimportant ones will probably be
% better, because it doesn't affect the overall power as much. 

% OPTION #1: Zeroing
% % Try shrinking at 100 different amounts, and recomputing the effects
% for shrink_factor = 1:10
% end

% Normalize effects to be more or less 1

fprintf('Picking magic ad-hoc power...\n');

% Trade a drop in performance for as much sparsity as we can handle
for aa = 0.00:0.01:5    
    % OPTION #2: Scale all coefficients
    for ii = 1:length(mod_idxs)
        idx = mod_idxs{ii};
        for jj = 1:length(STACK{idx})
            for kk = 1:length(coefs{ii}{jj}(:))
                STACK{idx}{jj}.coefs(kk) = 100.^aa * effects{ii}{jj}(kk).^aa * STACK{idx}{jj}.coefs(kk);
            end
        end
    end
    
    calc_xxx(idx);
    perf = META.perf_metric();
    fprintf('%f: Scored: %f\n', aa, perf);
    
    if (abs(init_perf - perf)/init_perf > N_fraction)
        fprintf('Stopping, because performance difference reached.\n');
        break
    end
    STACK = stack_cache;
end
