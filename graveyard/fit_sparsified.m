function [termcond, n_iters] = fit_sparsified(fitter, n_jacks, n_sparsity_loops, shrinkstyle)
% [termcond, n_iters] = fit_sparsified(fitter, n_jacks, n_sparsity_loops, shrinkstyle)
%
% Attempts to find parameters are optimally sparse.
%
% REQUIREMENTS:
% 1. Assumes that global variable XXX{2}.dat exists and has respfiles
%    loaded.
% 
% ALGORITHM:
% An outer loop or function tries different levels of sparsity, trying more
% or less sparsity as is appropriate.
%
% For each level of sparseness weight n_jacks jackknifes of the data are
% created, used to train n_jacks models, and the predictions of these
% models on the held-out data are then recorded. This creates a
% "prediction" vector exactly the same length as the original data.
%
% This process terminates after n_innerloops has expired. 10-20 is a good
% default range and will try sparsity precisions to roughly half an order
% of magnitude.

global STACK XXX META;

n_iters = 0;

if ~exist('fitter', 'var')
    fitter = @() fit_boost(25);
end

if ~exist('n_jacks', 'var')
    n_jacks = 5;
end

if ~exist('n_sparsity_loops', 'var')
    n_sparsity_loops = 10;
end

if ~exist('shrinkstyle', 'var')
    shrinkstyle = 'mean'; % Options: 'mean', 'james', 'stephen', 'besttrain', or 'bestholdout'
end
phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;  
    return 
end

[fit_start_depth, fit_end_depth] = find_fit_start_depth(STACK);

% Make sure the model structure is appropriate
% The loader modules must be OUTSIDE the fittable range, and 1 or more must exist
% The perf metric module must be OUTSIDE the fittable range, and exactly 1 must exist
areloaders = find(cellfun(@(x) isfield(x{1}, 'is_data_loader'), STACK));
areperfmetrics = find(cellfun(@(x) isfield(x{1}, 'is_perf_metric'), STACK));
if ~all(areloaders < fit_start_depth) || isempty(areloaders)
    error('Model structure not jackknife-able: data loaders must come before all modules with fittable fields.');
end
if ~all(areperfmetrics > fit_end_depth) || isempty(areperfmetrics)
    error('Model structure not jackknife-able: performance metric must come after all modules with fittable fields.');
end

% We will jackknife the data at the point where the fittables begin. 
% Therefore, we should cache the original value at that point.
cached_xxx = XXX;
cached_stack = STACK;

% Make sure all elements are the same length
for iii = 1:length(XXX{fit_start_depth}.training_set),
    sff = XXX{fit_start_depth}.training_set{iii};
    signals = fieldnames(XXX{fit_start_depth}.dat.(sff));  
    lengths = zeros(1, length(signals));
    for jjj = 1:length(signals)
        lengths = size(XXX{fit_start_depth}.dat.(sff).(signals{jjj}), 1);
    end
    if ~all(lengths == lengths(1))
        error('I cannot jackknife dat signals unless they are all exactly the same length.');
    end
end 

% ---------------------------------------------------------------------
% INNER LOOP FUNCTION
function [score, phi_jacked] = calc_jackknifed_prediction_score()
    xxx_jack = cell(1, n_jacks);
    stack_jack = cell(1, n_jacks);
    
    % Compute all the jackknifes, making copies of XXX and STACK
    for jj = 1:n_jacks
        fprintf('\nJackknife [%d/%d]\n', jj, n_jacks);
        
        % Start with the initial, full data set
        XXX{fit_start_depth}.dat = cached_xxx{fit_start_depth}.dat;        

        % Clear the test_set, since we are overwriting it temporarily
        XXX{fit_start_depth}.test_set = {};
        
        for ii = 1:length(XXX{fit_start_depth}.training_set)
            sf = XXX{fit_start_depth}.training_set{ii};
            
            nsf = [sf '_holdout']; % Fake training set stimfile            
            XXX{fit_start_depth}.test_set{ii} = nsf;
            
            sigs = fieldnames(XXX{fit_start_depth}.dat.(sf));
            for ss = 1:length(sigs),                                
                % Copy the training set data to the fake test set
                XXX{fit_start_depth}.dat.(nsf).(sigs{ss}) = XXX{fit_start_depth}.dat.(sf).(sigs{ss});
    
                % FIXME: This IF statement is a hack that is not generalizable!!!
                % I am NANing out everything in the respavg but not in
                % other signals because doing just the respavg avoids creating gaps in the
                % prediction caused by FIR filter passing over NANs
                if strcmp(sigs{ss}, 'respavg')                    
                    datdim=size(XXX{fit_start_depth}.dat.(sf).(sigs{ss}));
                    mask=ones(datdim(1),datdim(2));
                    jackidx=floor(linspace(1,1+datdim(1)*datdim(2),n_jacks+1));
                    mask(jackidx(jj):jackidx(jj+1)-1)=0;
                    mask=repmat(mask,[1 1 datdim(3:end)]);
                
                    % NAN out the held-out data in the training set stimfile
                    XXX{fit_start_depth}.dat.(sf).(sigs{ss})(mask==0) = NaN;
                
                    % Nan out everything but the held-out data in the test set
                    XXX{fit_start_depth}.dat.(nsf).(sigs{ss})(mask==1) = NaN;

                    % Old Way
                    %jackidx = floor(linspace(1, 1+size(XXX{fit_start_depth}.dat.(sf).(sigs{ss}), 1), n_jacks+1))
	                %XXX{fit_start_depth}.dat.(sf).(sigs{ss})(jackidx(jj):jackidx(jj+1)-1,:,:) = NaN;
                    %mask = ones(size(XXX{fit_start_depth}.dat.(sf).(sigs{ss})));
                    %mask(jackidx(jj):jackidx(jj+1)-1,:,:) = 0;
                    %XXX{fit_start_depth}.dat.(nsf).(sigs{ss})(mask==1) = NaN;
                end
            end
        end
        
        % Fit that jackknife
        unpack_fittables(phi_init);
        calc_xxx(fit_start_depth); 
        fitter();
                
        % Correct the sign of the model
        verify_model_polarity();
        
        % Save the STACK and XXX into the jackknife data structure
        xxx_jack{jj} = XXX;
        stack_jack{jj} = STACK; 
    end
            
%     % Plot some debugging crap
%     sel = [];    
%     sel.chan_idx = 1;
%     sel.stim_idx = 1;
     for ii = 1:n_jacks
%         %sel.stimfile = xxx_jack{ii}{5}.training_set{ii};
%         %figure; stack_jack{ii}{4}{1}.plot_fns{1}.fn(sel, stack_jack{ii}(1:4), xxx_jack{ii}(1:5));
%         %figure; stack_jack{ii}{4}{1}.plot_fns{1}.fn(sel, stack_jack{ii}(1:4), xxx_jack{ii}(1:5));
%         
         STACK = stack_jack{ii};
         XXX = xxx_jack{ii};
         plot_model_summary();       
%         %narf_modelpane
     end     
          
    % Compute which jackknifes had the best training and holdout scores
    bst_train_score = [];
    bst_train_idx = 0;
    bst_holdout_score = [];
    bst_holdout_idx = 0;
    for ii = 1:n_jacks
        STACK = stack_jack{ii};
        XXX = xxx_jack{ii};        
        [score, ~] = META.perf_metric();
        if isempty(bst_train_score) || bst_train_score > score
            bst_train_score = score;
            bst_train_idx = ii;
        end
        
        if isempty(bst_holdout_score) || (isfield(XXX{end}, 'score_test_corr') && ...
                                        bst_holdout_score > XXX{end}.score_test_corr)
            bst_holdout_score = XXX{end}.score_test_corr; % FIXME: Unfortunate Special case!
            bst_holdout_idx = ii;
        end
    end
        
    % Reload the original stack and xxx
    XXX = cached_xxx;
    STACK = cached_stack;
    
    % Merge the output of the jackknifed data sets at the input to the
    % performance metric module    
    idx = areperfmetrics(1);
        
    % Replace the training set with the held-out data
    for ii = 1:length(XXX{idx}.training_set)
        sf = XXX{idx}.training_set{ii};
        nsf = [sf '_holdout'];
        sigs = fieldnames(XXX{idx}.dat.(sf));
        for ss = 1:length(sigs),
            jackidx = floor(linspace(1, 1+size(XXX{idx}.dat.(sf).(sigs{ss}), 1), n_jacks+1));
            for jj = 1:n_jacks,
                % Beware all ye who look into the eye of madness:
                XXX{idx}.dat.(sf).(sigs{ss})(jackidx(jj):jackidx(jj+1)-1,:,:) = xxx_jack{jj}{idx}.dat.(nsf).(sigs{ss})(jackidx(jj):jackidx(jj+1)-1,:,:);
            end
        end
    end
    
    % Compute the performance from here to the end
    calc_xxx(idx);    
    [score, ~] = META.perf_metric();
    
    phi_jacks = cell2mat(cellfun(@pack_fittables, stack_jack, ...
                                 'UniformOutput', false));
    m = length(phi_init);    
    mu = mean(phi_jacks, 2);    
    sigma_sq = var(phi_jacks, [], 2);
        
    if strcmp('besttrain', shrinkstyle)
        phi_jacked = phi_jacks(:, bst_train_idx);
    elseif strcmp('bestholdout', shrinkstyle)
        phi_jacked = phi_jacks(:, bst_holdout_idx);
    elseif strcmp('mean', shrinkstyle)
        phi_jacked = mu;
    elseif strcmp('stephen', shrinkstyle)
        % Use Stephen's shrinkage equation
        phi_jacked = real(mu .* sqrt(1 - ((sqrt(sigma_sq) / sqrt(n_jacks)) ./ mu).^2));
        phi_jacked(isnan(phi_jacked)) = 0;
    elseif strcmp('james', shrinkstyle) 
        if m < 3
            fprintf('WARNING: A James Stein Estimator is only better if there are 3 or more params. Using simple mean.');
            phi_jacked = mu;
        else
            prior = phi_init;
            phi_jacked = prior + (mu - prior) .* (  1 - (((m-2) * (sigma_sq / n_jacks)) / sum((mu - prior).^2))); 
        end
    else
        error('shrinkstyle can only be ''besttrain'', ''mean'', ''stephen'', or ''james''.');
    end        
end

% ---------------------------------------------------------------------
% OUTER LOOP FUNCTIONS

% Sparsity search constants
sparsity_weight_max = 10^2;
sparsity_weight_min = 10^-8;

function [best_sparsity, best_phi] = linear_search ()
    sparsity_weights = logspace(log10(sparsity_weight_min), ...
                                log10(sparsity_weight_max), ... 
                                n_sparsity_loops);
	
    % Compute non-sparsity weighted starting point
    [best_score, ~] = META.perf_metric();
    best_phi = phi_init;
    best_sparsity = 0.0;
    
    for ii = 1:length(sparsity_weights)        
        META.sparsity_weight = sparsity_weights(ii);
        [score, phi] = calc_jackknifed_prediction_score();
        if score < best_score
            best_score = score;
            best_phi = phi;
            best_sparsity = sparsity_weights(ii);
        end
        
        fprintf('Sparsity level %e had score %f\n', ... 
            META.sparsity_weight, score);
    end
end

function [best_sparsity, best_phi] = bisection_search()
    % TODO: Exploit a better search algorithm than linear!
end

% Start one of the outer loops
[best_sparsity, best_phi] = linear_search();
%% [best_sparsity, best_phi] = bisection_search();

% Set the optimal sparseness weight level and best parameters found
META.sparsity_weight = best_sparsity;
unpack_fittables(best_phi);
calc_xxx(fit_start_depth); 

end