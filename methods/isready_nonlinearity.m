function isready = isready_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    if all(isfield(x, {'dat'}))
        isready = true;
        for sf = fieldnames(x.dat)', sf=sf{1};
            isready = isready && ...
                      all(isfield(x.dat.(sf), {'pp_stim', ...
                                               'raw_stim',...
                                               'raw_stim_time', ...
                                               'raw_stim_fs'}));
        end
    else
        isready = false;
    end
end
