function cache_summaries(cellid, skip_existing)
% Scans 'NARF_SAVED_MODELS_PATH/cellid/*.mat' for any model files.
% Loads the models, rebuilds the XXX path all the way through, and 
% caches quantities which are applicable to analysis in a 'model summary'
% file, which is stored as 'NARF_SAVED_ANALYSIS_PATH/cellid_summary.mat'.
%
% This is basically a way to build a cache so that later analyses can be
% performed without having to slowly load and process each model file. 
%
% ARGUMENTS:
%     CELLID          Self explanatory
%     SKIP_EXISTING   When true, filenames already in the analysis file are
%                     not loaded, allowing a cache to be rebuilt very 
%                     quickly to only include new models not yet cached.
%                     Defaults to false.

global NARF_SAVED_MODELS_PATH NARF_SAVED_ANALYSIS_PATH;

if nargin < 2
    skip_existing = false;
end

analysis_file = [NARF_SAVED_ANALYSIS_PATH filesep cellid '_summary.mat'];

if ~skip_existing && exist(analysis_file, 'file') == 2
    fprintf('Loading existing analysis file.\n');
    summary = getfield(load(analysis_file, 'summary'), 'summary');
else
    summary = {};
end

modelfiles = dir2cell([NARF_SAVED_MODELS_PATH filesep cellid filesep '*.mat']);
modelpaths = cellfun(@(x) [NARF_SAVED_MODELS_PATH filesep cellid filesep x], ...
                      modelfiles, 'UniformOutput', false);
                  
for ii = 1:length(modelfiles)
    mf = modelfiles{ii};
    mp = modelpaths{ii};
    fprintf('Loading model file %s\n', mf);
    load_model_stack(mp);
    summary{ii} = summarize_model();
end

save(analysis_file, 'summary');