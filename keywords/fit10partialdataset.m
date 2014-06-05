function fit10partialdataset()

global STACK XXX;

nmse();

    function fn = make_subfitter(del)
        function [a,b,c,d] = subfitter(prev_opts)
            
            % Detect whether the fittables are in a FIR block or not
            module_being_fit = '';
            for kk = 1:length(STACK)
                if isfield(STACK{kk}{1}, 'fit_fields') && ~isempty(STACK{kk}{1}.fit_fields)
                    module_being_fit = STACK{kk}{1}.name;
                    break;
                end
            end
            
            if strcmp(module_being_fit, 'fir_filter') || ...
                    strcmp(module_being_fit, 'weight_channels')
                if exist('prev_opts', 'var')
                    [a,b,c,d] = fit_boo(prev_opts);
                else
                    [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', del, ...
                        'StopAtStepNumber', 1, ...
                        'StepGrowth', 1.3);
                end
            else
                if exist('prev_opts', 'var')
                    [a,b,c,d] = fit_scaat(prev_opts);
                else
                    [a,b,c,d] = fit_scaat('StopAtAbsScoreDelta', del, ...
                        'StopAtStepNumber', 1);
                end
            end
        end
        
        fn = @subfitter;
        
    end

% Initialization: If FIR filters are all zero, initialize them randomly
[~, mod_idxs] = find_modules(STACK, 'fir_filter');
for ii = 1:length(mod_idxs)
    for jj = 1:length(STACK{mod_idxs{ii}})
        if isfield(STACK{mod_idxs{ii}}{jj}, 'fit_fields') && ...
                any(strcmp('coefs', STACK{mod_idxs{ii}}{jj}.fit_fields)) && ...
                all(0 == STACK{mod_idxs{ii}}{jj}.coefs(:))
            STACK{mod_idxs{ii}}{jj}.coefs = normrnd(0, 10^-3, size(STACK{mod_idxs{ii}}{jj}.coefs));
            fprintf('=====> Initializing FIR coefs to random numbers!\n');
        end
    end
end


% Split the training dataset in two parts.
r = size(XXX{2}.dat.(cell2mat(XXX{1}.training_set)).respavg,2);
inds = reshape(1:r,10,r/10);
mid = round(size(inds,1)/2);
testing_inds = mat2vec(inds(mid,:));
training_inds = mat2vec(inds([1:(mid-1) (mid+1):10],:));
training_XXX = XXX;
testing_XXX = XXX;
training_XXX{2}.dat.(cell2mat(XXX{1}.training_set)).respavg = ...
    XXX{2}.dat.(cell2mat(XXX{1}.training_set)).respavg(:,training_inds);
training_XXX{2}.dat.(cell2mat(XXX{1}.training_set)).stim = ...
    XXX{2}.dat.(cell2mat(XXX{1}.training_set)).stim(:,training_inds, :);
training_XXX{2}.dat.(cell2mat(XXX{1}.training_set)).resp = ...
    XXX{2}.dat.(cell2mat(XXX{1}.training_set)).resp(:,training_inds, :);
testing_XXX{2}.dat.(cell2mat(XXX{1}.training_set)).respavg = ...
    XXX{2}.dat.(cell2mat(XXX{1}.training_set)).respavg(:,testing_inds);
testing_XXX{2}.dat.(cell2mat(XXX{1}.training_set)).stim = ...
    XXX{2}.dat.(cell2mat(XXX{1}.training_set)).stim(:,testing_inds, :);
testing_XXX{2}.dat.(cell2mat(XXX{1}.training_set)).resp = ...
    XXX{2}.dat.(cell2mat(XXX{1}.training_set)).resp(:,testing_inds, :);

% one first run for nothing, to make sure the ...._XXX are well initialized
XXX = training_XXX;
calc_xxx(2);
training_XXX = XXX;
XXX = testing_XXX;
calc_xxx(2);
testing_XXX = XXX;


% now, a first rough search
XXX = training_XXX;
fit_boo('StopAtAbsScoreDelta', 10^-2, 'StepGrowth', 1.3);
training_XXX = XXX;

% compute the score in the testing sub-dataset
XXX = testing_XXX;
calc_xxx(2);
testing_XXX = XXX;
s = testing_XXX{end}.score_train_nmse;

% Now gradually shrink the stopping criterion, stopping when checking the
% (inner) testing dataset reveals overfitting for at least n consecutive
% set of iterations

scale=10^1;
% stop_at=10^-6;
stop_at=10^-6;
% group_size = 2; % this is the size of a "set of iterations". Checking with
% the testing dataset occurs only _after_ optimizing for
% this number of iterations.
decrease_consecutive_max = 2; % that's the maximal number of times that we
% are allowed to decrease
decrease_consecutive_done = 0;

last_good_stack = STACK;

while(scale > stop_at)
    %     n_iters = 0;
    
    %     while (n_iters<group_size)
    %     while (1)
    
    XXX = training_XXX;
    [term_cond, term_score, n_iters] = fit_iteratively(make_subfitter(scale), ...
        create_term_fn('StopAtAbsScoreDelta', scale)) %, ...
    %             'StopAtStepNumber', group_size));
    training_XXX = XXX;

    previous_s = s;
    % TODO: compute the score in the testing dataset, and put it in the
    % variable 's'
    
    % compute the score in the testing sub-dataset
    XXX = testing_XXX;
    calc_xxx(2);
    s = testing_XXX{end}.score_train_nmse;
    
    if previous_s < s
        % the testing dataset reveals overfitting => increase counter
        decrease_consecutive_done = decrease_consecutive_done + 1;
%         fprintf(1,'BOUH - this is getting worse\n');
%         fprintf(1,'     5555555555555555555555555555555555555555\n');
%         fprintf(1,'     5555555555555555555555555555555555555555\n');
%         fprintf(1,'     5555555555555555555555555555555555555555\n');
%         fprintf(1,'     55555555555555555555 testing NON-improvement = %d 555555555555555\n', previous_s - s);
%         fprintf(1,'     5555555555555555555555555555555555555555\n');
%         fprintf(1,'     5555555555555555555555555555555555555555\n');
%         fprintf(1,'     5555555555555555555555555555555555555555\n');
    else
        % last iteration was an improvement => reset counter and store
        % the STACK
        decrease_consecutive_done = 0;
        last_good_stack = STACK;
%         fprintf(1,'YAY\n');
%         fprintf(1,'     11111111111111111111111111111111\n');
%         fprintf(1,'     1111111111 testing improvement = %d 11111111111\n', previous_s - s);
%         fprintf(1,'     11111111111111111111111111111111\n');
        
    end
    
    if decrease_consecutive_done == decrease_consecutive_max
        % we went through "decrease_consecutive_max" bad steps
        % => put back the latest good stack, and then exit this madness
%         STACK = last_good_stack;
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
%         fprintf(1,'     99999999999999999999999999999999999999999999999999999\n');
        
        fprintf(1,'I should have quit!\n');
%         return
    end
    
    %     end
    
    
    scale = scale * 0.7; % Very conservative: 0.8. Probably 0.5 or even 0.25 is fine.
%     scale = scale * 0.1; % Very conservative: 0.8. Probably 0.5 or even 0.25 is fine.
end

end
