function [termcond, n_iters] = fit_split_simply(fitter, splitter, unifier)
% [termcond, n_iters] = fit_split(fitter, splitter, unifier);

% Algorithm Overview
%   1. Does not do initialization. If you want to fit everything all
%      together, you should do that before even calling this function.
%   2. Does not actually use splitters in each module during training.
%      Instead, this algorithm just fits each split case independently,
%      then joins everything together at the end and adds the splitter and
%      unifier to each module if it has a 

global STACK XXX META;

n_iters = 0;

if ~exist('fitter', 'var')
    fitter = @fit_boost;
end

if ~exist('splitter', 'var')
    splitter = @split_by_respfile;
end

if ~exist('unifier', 'var')   
    unifier = @unify_respfiles;
end

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;  
    return
end

[fit_start_depth, fit_end_depth] = find_fit_start_depth(STACK);

% Make sure the model structure is appropriate
% 1. The loader modules must be OUTSIDE the fittable range, and 1 or more must exist
% 2. The perf metric module must be OUTSIDE the fittable range, and exactly 1 must exist
% 3. No existing splitter/unifiers must be in the model
areloaders = find(cellfun(@(x) isfield(x{1}, 'is_data_loader'), STACK));
areperfmetrics = find(cellfun(@(x) isfield(x{1}, 'is_perf_metric'), STACK));
has_splitters = find(cellfun(@(x) isfield(x{1}, 'splitter'), STACK), 1);

if ~all(areloaders < fit_start_depth) || isempty(areloaders)
    error('Model structure not splittable: data loaders must come before all modules with fittable fields.');
end
if ~all(areperfmetrics > fit_end_depth) || isempty(areperfmetrics)
    error('Model structure not splittable: performance metric must come after all modules with fittable fields.');
end
if ~isempty(has_splitters)
    error('fit_split_simply() relies on all modules NOT being splittable.');
end

cached_STACK = STACK;
cached_XXX = XXX;

% Split the XXX data structure according
split_XXXs = splitter(XXX(1:fit_start_depth)); 
n_splits = length(split_XXXs);

% Make split versions of the stack
split_STACKs = {};
for ii = 1 :n_splits
    split_STACKs{ii} = STACK;
end

% Train on each of the splits
for ii = 1:n_splits
    fprintf('\nFitting split [%d/%d]\n', ii, n_splits);
    
    % Reset the STACK and XXX
    XXX = split_XXXs{ii};
    STACK = split_STACKs{ii}; 
    
    calc_xxx(fit_start_depth);     
    fitter();
    
    split_XXXs{ii} = XXX;
    split_STACKs{ii} = STACK;
end

STACK = cached_STACK;
XXX = cached_XXX;

% Merge all the splits together
for ii = fit_start_depth:length(STACK)
    %if isfield(STACK{ii}{1}, 'fit_fields')
        mm = STACK{ii}{1};
        if isfield(mm, 'is_splittable') && mm.is_splittable
            for jj = 1:n_splits
                STACK{ii}{jj} = split_STACKs{jj}{ii}{1};            
                STACK{ii}{jj}.splitter = splitter;
                STACK{ii}{jj}.unifier = unifier;
            end
        end
        %end
end

% Recalc once to create the updated performance metric
calc_xxx(fit_start_depth); 

end