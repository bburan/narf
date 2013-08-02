function [stacks, metas, x0s] = load_model_batch(batch, cellids, modelnames, ...
                                stack_extract_fn, meta_extract_fn, x0_extract_fn)
% [stacks, metas, x0s] = load_model_batch(batch, cellids, modelnames)
%
% Loads every NARF model in BATCH with a matching cellid and modelanme 
% found in CELLIDS and MODELNAMES. This function is intended to be used
% when extracting a bunch of parameters from completed NARF simulations.
%
% If you don't want to extract the entire stack, meta, or x0 matrix for 
% each model, you may provide three additional extra arguments, which are
% functions that will be applied to each stack to extract the desired 
% parameters that will then be bundled into 
%
% ARGUMENTS
%    filepath    The absolute path to where the NARF .mat file is found
% 
% RETURNS:
% Three 2x2 cell arrays are returned, with dimension 1 being the modelnames
% in the order presented in MODELNAMES, and dimension two being the cellids
% as they appear in CELLIDs. 
%
%    STACKS     A 2x2 cell array of the STACKs of every model. 
%    METAS      A 2x2 cell array of the METAs of every model. 
%    X0S        A 2x2 cell array of the XXX{1}'s of every model.

global STACK XXX META;

if ~exist('stack_extract_fn', 'var')
    stack_extract_fn = @(x) x;
end
if ~exist('meta_extract_fn', 'var')
    meta_extract_fn = @(x) x;
end
if ~exist('x0_extract_fn', 'var')
    x0_extract_fn = @(x) x;
end

n_models = length(modelnames);
n_cellids = length(cellids);
stacks = cell(n_models, n_cellids);
metas = cell(n_models, n_cellids);
x0s = cell(n_models, n_cellids);

dbopen;
for mi = 1:n_models
    model = modelnames{mi};
    fprintf('Loading %d modelfiles matching: %s\n', n_cellids, model);
    for ci = 1:n_cellids
        cellid = cellids{ci};
        
        sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(batch) ''];
        sql = [sql ' AND cellid="' cellid '"'];
        sql = [sql ' AND modelname="' model '"'];
        results = mysql(sql);  
        
        if isempty(results)
            stacks{mi, ci} = {};
            metas{mi, ci} = {};
            x0s{mi, ci} = {};
            continue;
        elseif (length(results) > 1)
            error('Multiple DB Results found: for %s/%s\n', cellid, model);
        end
        
        load_model(char(results(1).modelpath));
        stacks{mi, ci} = stack_extract_fn(STACK);
        metas{mi, ci} = meta_extract_fn(META);
        x0s{mi, ci} = x0_extract_fn(XXX{1});
    end
end