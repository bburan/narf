function fh = plot_best_models(cellid, summaries, n, savetodisk, sortfield)
% Uses compare_models() to view the best-performing models. Assumes that a
% cellid summary file has been built already. 
% 
% ARGUMENTS:
%    CELLID     Self-explanatory
%    N          Number of models to display.
%               Defaults to 5 if not provided.
%    SAVETODISK If true, a PNG file will be created in the analysis dir.
%               Defaults to false if not provided.
%    SORTFIELD  The field by which to sort the models.
%               Defaults to 'score_test_corr' if not provided.
%
% Returns a handle to the window that just opened. If savetodisk is true,
% the window will automatically be closed and the returned value is NaN. 

global NARF_SAVED_ANALYSIS_PATH;

if nargin < 3
    n = 5;
end

if nargin < 4
    savetodisk = false;
end

if nargin < 5
    sortfield = 'score_test_corr';
end

sr = sort_by_field(summaries, sortfield);

best = sr(max(1, end-n+1):end);

filepaths = extract_field(best, 'modelpath');

fh = compare_models(filepaths);

savethefig(fh, cellid, sprintf('%s_best%d.png', cellid, n));

end