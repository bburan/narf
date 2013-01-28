function fh = plot_best_models(cellid, n, savetodisk, sortfield)
% Uses compare_models() to view the best-performing models. 
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

if nargin < 2
    n = 5;
end

if nargin < 3
    savetodisk = false;
end

if nargin < 4
    sortfield = 'score_test_corr';
end
   
summary_file = [NARF_SAVED_ANALYSIS_PATH filesep cellid '_summary.mat'];

if exist(summary_file, 'file') ~= 2
    fprintf('Summary file not found: %s\n', summary_file);
    return
end

summary = getfield(load(summary_file, 'summary'), 'summary');

% Strip out any empty summary
smry = {};
for ii = 1:length(summary)
    if ~isempty(summary{ii})
        smry{end+1} = summary{ii};
    end
end

sr = sort_by_field(smry, sortfield);

best = sr(max(1, end-n+1):end);

filenames = extract_field(best, 'modelname');

fh = compare_models(filenames);
pngfile = [NARF_SAVED_ANALYSIS_PATH filesep cellid '.png'];
if savetodisk
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'InvertHardcopy','off');
    print(fh, pngfile, '-dpng');    
    close(fh);
    fh = nan;
end

end