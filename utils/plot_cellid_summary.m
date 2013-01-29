function fhs = plot_cellid_summary(cellid, summaries, savetodisk)
% Creates several figures which describe, and saves them to disk if
% SAVETODISK is true. Assumes that a cellid summary file has been built.
% 
% ARGUMENTS:
%    CELLID     Self-explanatory
%    SUMMARIES  Cell array of summary structs.
%    SAVETODISK If true, a PNG file will be created in the analysis dir.
%               Defaults to false if not provided.
%
% Returns cell array of handles to the windows that were opened. 
% If savetodisk is true, these windows will automatically be closed and
% the returned value is NaN. 

if nargin < 3
    savetodisk = false;
end

fhs = {};

function append_or_save(fh, filename)
    if savetodisk
        % Save to disk and close
        savethefig(fh, filename);
    else 
        fhs{end+1} = fh;
    end
end

% ------------------------------------------------------------------------
% TOP N MODEL SUMMARIES
n = 5;
sr = sort_by_field(summaries, 'score_test_corr');
best = sr(max(1, end-n+1):end);
filepaths = extract_field(best, 'modelpath');
fh = compare_models(filepaths);
append_or_save(fh, sprintf('%s_best%d', cellid, n));

% ------------------------------------------------------------------------
% TOKEN SCATTER PLOTS

% Build a list of summary tokens and a structure to store them in
tokens = cell(1, 10); % FIXME
for ii = 1:length(summaries)
    s = summaries{ii};
    toks = tokenize_modelname(s.modelname);
    
    for jj = 1:length(toks)
        tokens{jj} = cat(2, tokens{jj}, toks{jj});
    end  
end

% Condense tokens
for ii = 1:length(toks)
    tokens{ii} = unique(tokens{ii});
end

% Now generate scatter and bar plots for each of those
for ii = 2:4
    fh = plot_token_performance(summaries, [cellid sprintf(', token group %d', ii)], ii, tokens{ii});
    append_or_save(fh, sprintf('%s_tok%d', cellid, ii));
end

% ------------------------------------------------------------------------

end