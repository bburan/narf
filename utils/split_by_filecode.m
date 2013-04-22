function xxxs = split_by_filecode(xxx)
% {xxx1, xxx2, etc} = splitter(XXX)
%
% Returns a cell array of XXX structs with only one respfile per final
% cell.

estfiles = xxx{end}.training_set;
valfiles = xxx{end}.test_set;

xxxs = cell(1, length(estfiles));

% If there are no filecodes, put everything together
if isfield(xxx{1}, 'filecode')
    filecodes = xxx{1}.filecode;
else
    xxxs = {xxx};
    return
end

% Cut off everything after the first character of the filecode
for ii = 1:length(filecodes)
    tmp = filecodes{ii};
    tmp = tmp(1); % Just the first character
    filecodes{ii} = tmp;
end

% Otherwise, group by filecode
unique_codes = unique(filecodes);
for ii = 1:length(unique_codes);
    fc = unique_codes{ii};
    
    matches = strcmp(fc, filecodes);   
    
    efs = estfiles(matches);
    vfs = valfiles(matches); 
    
    tmp = xxx;
    tmp.dat = [];
    tmp.training_set = {};
    tmp.test_set = {};
    
    for jj = 1:length(fs)
        f = fs{jj};
        tmp.dat.(esf) = xxx{end}.dat.(esf);
        tmp.dat.(vsf) = xxx{end}.dat.(vsf);
        tmp.training_set{end+1} = esf;
        tmp.test_set{end+1}     = vsf;
    end
    
    xxxs{ii}{end+1} = tmp;
end
