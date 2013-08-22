function xxxs = split_by_filecodeP(xxx)
% {xxx1, xxx2, etc} = splitter(XXX)
%
% Returns a cell array of XXX structs with only one respfile per final
% cell.

% If there are no filecodes, put everything together
if isfield(xxx{1}, 'filecodes')
    filecodes = xxx{1}.filecodes;
else
    xxxs = {xxx};
    return
end

% Cut off everything after the first character of the filecode
for ii = 1:length(filecodes)
    tmp = filecodes{ii};
    tmp = tmp(1); % Just the first character
    if tmp=='P',
        filecodes{ii} = tmp;
    end
end

unique_codes = unique(filecodes);
xxxs = cell(1, length(unique_codes));
estfiles = xxx{end}.training_set;
valfiles = xxx{end}.test_set;

% Otherwise, group by filecode
for ii = 1:length(unique_codes);
    fc = unique_codes{ii};
    
    matches = strcmp(fc, filecodes);   
    
    efs = estfiles(matches);
    vfs = valfiles(matches(matches<=length(valfiles))); 
    
    xxxs{ii} = xxx;
    xxxs{ii}{end}.dat = [];   
    xxxs{ii}{end}.training_set={};
    xxxs{ii}{end}.test_set={};
    xxxs{ii}{end}.filecodes={};
    
    for jj = 1:length(efs)
        f = efs{jj};
        xxxs{ii}{end}.dat.(f) = xxx{end}.dat.(f);
        xxxs{ii}{end}.training_set{jj} = f;
        xxxs{ii}{end}.filecodes{jj} = fc;
    end
    for jj = 1:length(vfs)
        f = vfs{jj};
        xxxs{ii}{end}.dat.(f) = xxx{end}.dat.(f);
        xxxs{ii}{end}.test_set{jj} = f;
    end
    
end

