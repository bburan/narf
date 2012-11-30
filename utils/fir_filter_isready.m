function isready = fir_filter_isready(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    if all(isfield(x, {'dat'}))
        isready = true;
        for sf = fieldnames(x.dat)', sf=sf{1};
            isready = isready && ~isfield(x.dat.(sf), 'lf_stim') && ...
                      all(isfield(x.dat.(sf), {'ds_stim', ...
                                               'ds_stim_time',...
                                               'ds_stim_fs'}));
        end
    else
        isready =false;
    end
end