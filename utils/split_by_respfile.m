function xxxs = split_by_respfile(xxx)
% {xxx1, xxx2, etc} = splitter(XXX)
%
% Returns a cell array of XXX structs with only one respfile per final
% cell.

estfiles = xxx{end}.training_set;
valfiles = xxx{end}.test_set;

%if length(estfiles) ~= length(valfiles)
%    error('I cannot split by respfiles unless there are an equal number of estimation and validation files.');
%end

xxxs = cell(1, length(estfiles));

for ii = 1:length(estfiles)
    esf = estfiles{ii};    
    xxxs{ii} = xxx;
    xxxs{ii}{end}.dat = [];   
    xxxs{ii}{end}.dat.(esf) = xxx{end}.dat.(esf);
    xxxs{ii}{end}.training_set = {esf};
    
    if length(valfiles) >= ii
        vsf = valfiles{ii};
        xxxs{ii}{end}.dat.(vsf) = xxx{end}.dat.(vsf);
        xxxs{ii}{end}.test_set     = {vsf};
    end
    
    % Filecodes are not guaranteed to exist, so we need to do this:
    % FIXME: In a later version of baphy/narf, we should make explicit what
    % will be provided and what not
    if length(estfiles) == length(xxx{end}.filecodes)
        xxxs{ii}{end}.filecodes    = xxx{end}.filecodes(ii);
    else
        xxxs{ii}{end}.filecodes    = xxx{end}.filecodes;
    end
end
