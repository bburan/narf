function xxxs = split_by_filecode(xxx)
% {xxx1, xxx2, etc} = splitter(XXX)
%
% Returns a cell array of XXX structs with only one respfile per final
% cell.

fns = fieldnames(xxx{end}.dat);

xxxs = cell(1, length(fns));

% If there are no filecodes, put everything together
if isfield(xxx{1}, 'filecode')
    filecodes = xxx{1}.filecode;
else
    xxxs = {xxx};
    return
end

% Otherwise, group by filecode
tmp = [];
unique_codes = unique(filecodes);
for ii = 1:length(unique_codes);
    fc = unique_codes{ii};
    
    matches = strcmp(fc, filecodes);   
    
    fs = fns(matches);

    tmp = xxx;
    tmp.dat = [];
    for jj = 1:length(fs)
        f = fs{jj};
        tmp.dat.(f) = xxx{end}.dat.(f);
    end
    xxxs{ii}{end+1} = tmp;

end
