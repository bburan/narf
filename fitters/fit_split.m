function [termcond, n_iters] = fit_split(fitter, splitter, unifier)
% [termcond, n_iters] = fit_segmented(fitter, splitter, unifier);

global STACK XXX;

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
has_splitters = find(cellfun(@(x) isfield(x{1}, 'splitter'), STACK));

if ~all(areloaders < fit_start_depth) || isempty(areloaders)
    error('Model structure not splittable: data loaders must come before all modules with fittable fields.');
end
if ~all(areperfmetrics > fit_end_depth) || isempty(areperfmetrics)
    error('Model structure not splittable: performance metric must come after all modules with fittable fields.');
end
if any(has_splitters)
    error('Model structure not splittable: No splitters/unifiers can exist before this.');
end

cached_xxx = XXX;
cached_stack = STACK;

% Split the input to the model after the loader
loader_idx = find(cellfun(@(x) isfield(x{1}, 'is_data_loader'), STACK));
xxx_loaded_splits = splitter(XXX(loader_idx+1));

n_splits = length(xxx_loaded_splits);
xxx_splits = cell(1, n_splits);
stack_splits = cell(1, n_splits);

for ii = 1:n_splits
    fprintf('\nSplit [%d/%d]\n', ii, n_splits);
    
    % Start with the initial, full data set
    XXX = xxx_loaded_splits{ii};
    STACK = cached_stack;
    
    calc_xxx(loader_idx+1); 
    
    fitter();
    
    xxx_splits{ii} = XXX;
    stack_splits{ii} = STACK;
end

% Reload the original stack and xxx
XXX = cached_xxx;
STACK = cached_stack;

% Merge all the splits together using the splitter/unifier paramset thing
for ii = fit_start_depth:length(STACK)
    for jj = 1:n_splits
        STACK{ii}{jj} = stack_splits{jj}{ii}{1};
        STACK{ii}{jj}.splitter = splitter;
        STACK{ii}{jj}.unifier = unifier;
    end
end

% Recalc once to create the updated performance metric
calc_xxx(fit_start_depth); 

end