function plot_batch_summary(batches)
% plot_batch_summary(batches)
%
% For each cellid/training set combination in each batch,
%    Makes a list of all matching model files
%    If there are no cached image files, it generates them
%    
% Finally, it makes
%    Generates heat maps of model performance for cellids

if nargin < 1
    batches = [240];
end

baphy_set_path;
narf_set_path; 

global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();

for bi = 1:length(batches)
    batchnum = floor(batches(bi));
    
    analysis_prefix = sprintf('a%d', batchnum);
    
    cells = request_celldb_batch(batchnum);
    
    [~, modelnames] = module_combinations(mm);
    
    for ii = 1:length(cells)        
        % Technically, fit_models already built a perfect cache, even if interrupted,
        % but rebuilding it can give us a little peace of mind so let's do it.
        summarize_cellid(cells{ii}.cellid, true);
        
        % Otherwise, load the model summaries for this cell from the cache
        sf = [NARF_SAVED_ANALYSIS_PATH filesep cells{ii}.cellid '_summary.mat'];
        summaries = load_summaries({sf});
        
        % Limit the summaries to only those model combinations found in mm
        summaries = only_named_summaries(summaries, modelnames);
        
        % Generate PNGs showing the best models, tokens for each cellid
        plot_cellid_summary(cells{ii}.cellid, summaries, true, analysis_prefix);
        
    end
    
    % Finally, display and save heat maps of test and training set performance
    summary_files = cellfun(@(x) [NARF_SAVED_ANALYSIS_PATH filesep x.cellid '_summary.mat'], ...
        cells, 'UniformOutput', false);
    summaries = load_summaries(summary_files);
    summaries = only_named_summaries(summaries, modelnames);
    
    plot_summaries(summaries, analysis_prefix, true);
end


end