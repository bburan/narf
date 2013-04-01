function S = load_summaries(summary_files)
% Concatenate all analysis summaries into a big cell array, and return it.

S = {};
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
        S{end+1} = s;
    end
end
