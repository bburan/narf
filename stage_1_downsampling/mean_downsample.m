function ret = mean_downsample(args)
% Default Parameters
params = [];
params.raw_freq = 100000;
params.ds_freq = 200;
params.ds_dim_num = 2;  % The dimension on which to downsample
params.pretty_name = 'mean(abs(stim)), binned resp';
params.editable_fields = {'raw_freq', 'ds_freq'};

% Overwrite the defaults with the arguments iff they are in editable_fields
if nargin == 1
    fns = fieldnames(args);
    for idx = 1:length(fns);
        if any(ismember(params.editable_fields, fns{idx}))
            params.(fns{idx}) = args.(fns{idx});
        end
    end
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function [ds_stim, ds_resp] = do_mean_downsample(parm, stim, resp)
    scale = ceil(parm.raw_freq / parm.ds_freq);
    ds_stim=conv_fn(abs(stim), parm.ds_dim_num, @mean, scale, 0);
    ds_resp=conv_fn(resp, parm.ds_dim_num, @(x) (parm.ds_freq * mean(x)), scale, 0);
end

% ------------------------------------------------------------------------
% Put the instance methods in params struct
params.downsamp_fn = @do_mean_downsample;
ret = params;
end