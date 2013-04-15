function [xs, ys, mns] = do_plot_prep(xxs, stimfile, stim_idxs, chan_idxs, xfield, yfield )
    xs = {};
    ys = {};
    mns = {};
    
    for ii = 1:length(xxs)

        % Skip unless all fields are in existence
        if  ~isfield(xxs{ii}.dat, stimfile) || ...
            ~isfield(xxs{ii}.dat.(stimfile), xfield) || ...
            ~isfield(xxs{ii}.dat.(stimfile), yfield)
            continue;
        end
        
        mns{end+1} = ['PS' num2str(ii)];
        xs{end+1} = xxs{ii}.dat.(stimfile).(xfield)(:);
        ys{end+1} = squeeze(xxs{ii}.dat.(stimfile).(yfield)(:, stim_idxs, chan_idxs(:)));
        
    end
    
end