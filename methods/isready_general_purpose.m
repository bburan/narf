function isready = isready_general_purpose(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    if all(isfield(x, {'dat'}))
        isready = true;
        for sf = fieldnames(x.dat)', sf=sf{1};
            isready = isready && ...
                      all(isfield(x.dat.(sf), {'stim', ...
                                               'stim_time'}));
        end
    else
        isready = false;
    end
end
