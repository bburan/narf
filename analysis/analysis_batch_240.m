function analysis_batch_240(instance_num, total_instances)
% A script to analyze every cell in BAPHY's batch #240.
% It is expected you would call this from a bash shell script
% and that you want to run multiple instances of matlab simultaneously.
% EXAMPLE:
%    matlab -r "analysis_batch_240(1,4)" &
%    matlab -r "analysis_batch_240(2,4)" &
%    matlab -r "analysis_batch_240(3,4)" &
%    matlab -r "analysis_batch_240(4,4)" &

if nargin < 2
    instance_num = 1;
    total_instances = 1;
end

baphy_set_path;
narf_set_path; 

global MODULES NARF_SAVED_ANALYSIS_PATH;
MODULES = scan_directory_for_modules();

analysis_prefix = 'a240';

cells = request_celldb_batch(240);

% Define how groups of modules should be combined to make many models
mm = {}; 
mm{1} = module_groups('env100');
mm{2} = module_groups('log2');
mm{3} = module_groups('firb');
mm{4} = module_groups('senl', 'npnl');
mm{5} = module_groups('fminlsq', 'sb', 'boost');
mm{6} = module_groups('mse');

[~, modelnames] = module_combinations(mm);

for ii = 1:length(cells)
    if mod(ii + instance_num, total_instances) ~= 0
        continue
    end
    fit_models(mm, cells{ii}.cellid, ...
                   cells{ii}.training_set, ...
                   cells{ii}.test_set);
       
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
