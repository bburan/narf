function ret = meanlog_downsample(args)
% Default Parameters
params = [];
params.raw_freq = 100000;
params.ds_freq = 200;
params.ds_dim_num = 2;  % The dimension on which to downsample
params.pretty_name = 'Mean of sqrt(abs(x))';
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
function ds_x = do_meanlog_downsample(parm, x)
    ds_x=conv_fn(sqrt(abs(x)), parm.ds_dim_num, @mean, ceil(parm.raw_freq / parm.ds_freq), 0);
end

% ------------------------------------------------------------------------
% Put the instance methods in params struct
params.downsamp_fn = @do_meanlog_downsample;
ret = params;
end