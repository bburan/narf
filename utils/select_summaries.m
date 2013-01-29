function [M, xl, yl] = select_summary(summary_files, extractor_fn, ...
    selector_fn, xsort_fn, ysort_fn, xfield, yfield)
% Builds a summary matrix selecting only certain models using SELECT_FN,
% extracting some data using EXTRACT_FN, and sorting the X and Y axes using
% the XSORTER_FN and YSORTER_FN functions' return values. Only the first
% two arguments are required; other arguments have sane defaults.
%
% Incomplete values will be marked with NaN.
%
% ARGUMENTS:
%     SUMMARY_FILES    Which cellid summary files to use. A summary file
%                      may contain data for multiple cellids, but usually
%                      there will be just models for a single cellid.
%     SELECTOR_FN      Predicate function to determine if a model should be
%                      included in the summary matrix. SELECTOR_FN will be
%                      passed a model summary structure.
%     EXTRACTOR_FN     A function to extract a field from each model
%                      summary structure (which was produced by the
%                      summarize_model function).
%     XSORT_FN         A function which determines the sort order of the X
%                      axis labels. XSORTER_FN is not passed individual
%                      model summary structures, but is rather passed the
%                      entire column of sorted summary structures.
%     YSORT_FN         Same but for Y axis so rows are passed to YSORTER_FN
%
% RETURNS:
%     M    Data matrix full of the outputs of EXTRACT_FN.
%     XL   Cell array of x-axis labels. (Always CELLIDs)
%     YL   Cell array of y-axis labels. (Always Modelnames).
%
% EXAMPLE:
% Select only files with '_fir_' in their model names, extract test scores,
% and sort by average training score of that model.
% TODO!
% Select only files with 'log' in their model names, extract fit_time,
% and sort by fitter.
% TODO!

S = {};

if nargin < 2
    error('select_summaries() requires at least 2 arguments.');
end

if nargin < 3
    selector_fn = @(s) true; % Default is to select everything
end

if nargin < 4
    xsort_fn = @(ss) ss{1}.modelname;
end

if nargin < 5
    ysort_fn = @(ss) ss{1}.cellid;
end

if nargin < 5
    xfield = 'modelname';
end

if nargin < 5
    yfield = 'cellid';
end

% Put all analysis summaries into a big struct array
for ii = 1:length(summary_files)
    sf = summary_files{ii};
    
    % Skip if the summary file doesn't exist
    if exist(sf, 'file') ~= 2
        fprintf('Skipping nonexistant summary file: %s\n', sf);
        continue;
    end
    
    summary = getfield(load(sf, 'summary'), 'summary');
    
    % Concatenate those summaries onto S if they are selected
    for ii = 1:length(summary)
        s = summary{ii};
        if isempty(s)
            continue;
        end
        if selector_fn(s)
            S{end+1} = s;
        end
    end
end

% A reverse lookup table builder to find index number based on value
% Only works for cell arrays in which every element is unique.
    function idxstruct = build_index_lookup(cellarray)
        idxstruct = [];
        for jj = 1:length(cellarray)
            idxstruct.(cellarray{jj}) = jj;
        end
    end

% Chunk the data into X and Y bins
xchunks = [];
ychunks = [];
for ii = 1:length(S)
    xf = S{ii}.(xfield);
    xf = regexprep(xf, '-', '_');
    
    if isfield(xchunks, xf)
        xchunks.(xf){end+1} = S{ii};
    else
        xchunks.(xf) = {S{ii}};
    end
    
    yf = S{ii}.(yfield);
    yf = regexprep(yf, '-', '_');  
    if isfield(ychunks, yf)
        ychunks.(yf){end+1} = S{ii};
    else
        ychunks.(yf) = {S{ii}};
    end   
end

% Sort the bins based on xsort_fn and ysort_fn
xk = fieldnames(xchunks);
yk = fieldnames(ychunks);
xv = cellfun(@(k) xsort_fn(xchunks.(k)), xk, ...
                        'UniformOutput', false);
yv = cellfun(@(k) ysort_fn(ychunks.(k)), yk, ...
                        'UniformOutput', false);
if all(isnumeric(xv))
    xv = cell2mat(xv);
end
if all(isnumeric(yv))
    yv = cell2mat(yv);
end
                    
[xvl, ox] = sort(xv);
[yvl, oy] = sort(yv);    
    
    function v = safeget(x, f)
        if isfield(x, f)
            v = getfield(x, f);
        else
            v = nan;
        end
    end

% Extract the labels from each bin
xl = cellfun(@(k) safeget(xchunks.(k){1}, xfield), xk(ox), 'UniformOutput', false);
yl = cellfun(@(k) safeget(ychunks.(k){1}, yfield), yk(oy), 'UniformOutput', false);                    
xl = regexprep(xl, '-', '_');
yl = regexprep(yl, '-', '_');

% Build the X and Y lookup tables
xi_lookup = build_index_lookup(xl);
yi_lookup = build_index_lookup(yl);

% Now that we know the width and height, we can allocate a NaN matrix
M = NaN * zeros(length(yl), length(xl));

% Replace the NaNs with data, when we have it.
for ii = 1:length(S)
    xf = S{ii}.(xfield);
    xf = regexprep(xf, '-', '_');
    yf = S{ii}.(yfield);
    yf = regexprep(yf, '-', '_');  
    M(yi_lookup.(yf), xi_lookup.(xf)) = extractor_fn(S{ii});
end

end
