function ret = fir_filter(args)
% A 1D FIR Filter is created for each input filter dimension.
% The total number of filter coefficients is: num_coefs*num_filts

% Default Parameters
params = [];
params.num_coefs = 20;
params.num_filts = 2;
params.coefs = [0 5 4 3 1 -1 -2 -2 -1 0 0 0 0 0 0 0 0 0 0 0;
                1 1 1 1 1  1  1  1  1 1 1 1 1 1 1 1 1 1 1 1];
params.filt_dim_num = 2;
params.pretty_name = 'FIR Filters (1st Order Indep.)';
params.editable_fields = {'num_coefs', 'num_filts'};

% Overwrite the defaults with the arguments iff they are in editable_fields
if nargin == 1
    fns = fieldnames(args);
    for idx = 1:length(fns);
        if any(ismember(params.editable_fields, fns{idx}))
            params.(fns{idx}) = args.(fns{idx});
        end
    end
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if ~isequal([params.num_filts params.num_coefs], size(params.coefs))
    params.coefs = zeros(params.num_filts, params_num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function pred = do_fir_filtering(parm, x)
    % Apply the FIR filters to corresponding downsampled stimuli to get the model prediction
    % Since it is linear, the prediction is just the sum of the filters
    % We assume that there are no second order terms combining elements of both filters
    preds = [];
    for filt_idx = 1:parm.num_filts
        preds = cat(3, preds, filter(parm.coefs(filt_idx,:), 1, ...
                                       squeeze(x(:,:,filt_idx)), [], ...
                                       parm.filt_dim_num));
    end
    pred = squeeze(sum(preds, 3));
end

function plot_fir_coefs(parm)
    for filt_idx = 1:parm.num_filts
        stem([1:parm.num_coefs], parm.coefs(filt_idx,:), pickcolor(filt_idx));
    end
    %setAxisLabelCallback(gca, @(t) (t-1), 'X');
    axis tight;
end

% ------------------------------------------------------------------------
% Put the instance methods in params struct
params.model_fn = @do_fir_filtering;
params.plot_coefs_fn = @plot_fir_coefs;
ret = params;
end