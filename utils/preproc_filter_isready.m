function isready = preproc_filter_isready(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    % We are ready iff the necessary fields exist for every data file
    % in the .dat substructure
    
    % TODO: We also need to check that load_stim_resps_from_baphy exists
    if all(isfield(x, {'dat'})) 
        isready = true;
        for sf = fieldnames(x.dat)', sf=sf{1};
            isready = isready && ~isfield(x.dat.(sf), 'pp_stim') && ...
                      all(isfield(x.dat.(sf), {'raw_stim', 'raw_stim_time', 'raw_stim_fs'}));
        end     
    else
        isready =false;
    end
end
