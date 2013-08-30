function [termcond, n_iters] = fit_shrinkified3(fitter, n_jacks, n_shrink_loops, weighted_avg, use_corr_during_jackknifes)
% [termcond, n_iters] = fit_shrinkified(fitter, n_jacks, n_shrink_loops,
% weighted_avg, use_corr_during_jackknifes)

global STACK XXX META;

n_iters = 0;

if ~exist('fitter', 'var')
    fitter = @fit_boost;
end

if ~exist('n_jacks', 'var')
    n_jacks = 10;
end

if ~exist('n_shrink_loops', 'var')   
    n_shrink_loops = 100;
end

if ~exist('weighted_avg', 'var')
    weighted_avg = false;
end

if ~exist('use_corr_during_jackknifes', 'var')
    use_corr_during_jackknifes = false;
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
cached_perf_metric = META.perf_metric;

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
function [stack_jack, xxx_jack] = calc_jackknifes()
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
                  
    % Reload the original stack and xxx
    XXX = cached_xxx;
    STACK = cached_stack;    
end

function [holdout_score, phi_jacked] = shrink_and_predict(stack_jack, xxx_jack, shrink_amount)           
    phi_jacks = cell2mat(cellfun(@pack_fittables, stack_jack, ...
                                 'UniformOutput', false));
                                
    XXX = cached_xxx;
    STACK = cached_stack;    
    % Replace the training set with the held-out data
%     for ii = 1:length(XXX{fit_start_depth}.training_set)
%         sf = XXX{fit_start_depth}.training_set{ii};
%         nsf = [sf '_holdout'];
%         sigs = fieldnames(XXX{fit_start_depth}.dat.(sf));
%         for ss = 1:length(sigs),
%             jackidx = floor(linspace(1, 1+size(XXX{fit_start_depth}.dat.(sf).(sigs{ss}), 1), n_jacks+1));
%             for jj = 1:n_jacks,
%                 % Beware all ye who look into the eye of madness:
%                 XXX{fit_start_depth}.dat.(sf).(sigs{ss})(jackidx(jj):jackidx(jj+1)-1,:,:) = xxx_jack{jj}{fit_start_depth}.dat.(nsf).(sigs{ss})(jackidx(jj):jackidx(jj+1)-1,:,:);
%             end
%         end
%     end
    
    m = length(phi_init);
    mu = mean(phi_jacks, 2);
    sigma_sq = var(phi_jacks, [], 2);    
    phi_jacked = real(mu .* sqrt(1 - shrink_amount*((sqrt(sigma_sq) / sqrt(n_jacks)) ./ mu).^2));
    phi_jacked(isnan(phi_jacked)) = 0;
        
    unpack_fittables(phi_jacked);
    calc_xxx(fit_start_depth); 
    [holdout_score, ~] = META.perf_metric();
    
    fprintf('Holdout_score/shrink_amount: %e / %e\n', holdout_score, shrink_amount);
    
end

% ---------------------------------------------------------------------
% OUTER LOOP FUNCTIONS

% Sparsity search constants
shrink_max = 10^5;
shrink_min = 10^-5;

function [best_phi] = linear_search ()
    shrinkage = logspace(log10(shrink_min), ...
                         log10(shrink_max), ... 
                         n_shrink_loops);
	
	META.sparsity_weight = 0.0;
   
    % Calculate the Jackknifes, which is slow
	[stack_jack_orig, xxx_jack_orig] = calc_jackknifes();

    % Compute non-sparsity weighted starting point
    best_phi = phi_init;
    best_holdout_score = NaN;
    
    if use_corr_during_jackknifes
        META.perf_metric = @pm_corr;
    end
    % Try a bunch of shrinkage levels (which is fast) and find the best
    for ii = 1:length(shrinkage)        
        [holdout_score, shrunk_phi_jack_avg] = shrink_and_predict(stack_jack_orig, xxx_jack_orig, shrinkage(ii));
        if isnan(best_holdout_score) || holdout_score < best_holdout_score
            best_holdout_score = holdout_score;
            best_phi = shrunk_phi_jack_avg;
        end
    end
    
    % Replace the perf metric
    META.perf_metric = cached_perf_metric;
end

[best_phi] = linear_search();

% Set the optimal sparseness weight level and best parameters found
META.sparsity_weight = 0.0;
unpack_fittables(best_phi);
calc_xxx(fit_start_depth); 

end