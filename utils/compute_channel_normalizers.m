function [rms, dc_offset] = compute_channel_normalizers(x, signal, force_positive)
    tstim=[];
    fns = fieldnames(x.dat);
    for ii = 1:length(fns),
        sf = fns{ii};
        [T, S, C] = size(x.dat.(sf).(signal));
        tstim = cat(1, tstim, reshape(x.dat.(sf).(signal),T*S,C));
    end
    if force_positive
        mm = nanmin(tstim);
    else
        mm = nanmean(tstim);
    end
    rms=nanstd(tstim);    
    dc_offset = mm;
    
end