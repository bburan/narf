function xxxs = split_by_respfile(xxx)
% {xxx1, xxx2, etc} = splitter(XXX)
%
% Returns a cell array of XXX structs with only one respfile per final
% cell.

fns = fieldnames(xxx{end}.dat);

xxxs = cell(1, length(fns));

for ii = 1:length(fns)
    sf = fns{ii};
    xxxs{ii} = xxx;
    xxxs{ii}{end}.dat = [];   
    xxxs{ii}{end}.dat.(sf) = xxx(end).dat.(sf);
end
