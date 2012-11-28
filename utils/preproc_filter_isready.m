function isready = preproc_filter_isready(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    % We are ready iff the necessary fields exist for every data file
    % in the .dat substructure
    
    % TODO: We also need to check that load_stim_resps_from_baphy exists
    if all(isfield(x, {'dat'})) 
        sfs = fieldnames(x.dat);
        isready = true;
        for idx = 1:length(sfs)
            isready = isready && ...
                      all(isfield(x.dat.(sfs{idx}), {'raw_stim', 'raw_stim_time', 'raw_stim_fs'}));
        end     
    else
        isready =false;
    end
end
