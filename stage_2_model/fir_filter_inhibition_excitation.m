function ret = fir_filter_inhibition_excitation(args)
% A pair of 1D FIR Filters is created for each input filter dimension.
% The total number of filter coefficients is: num_coefs*num_filts*2
% One filter has positive effect on the output (Excitation)
% The other can only have a negative effect on the output (Inhibition)

% Default Parameters
params = [];
params.num_coefs = 20;
params.num_filts = 2;
params.excit_sigmoid = [1 2 3];
params.inhib_sigmoid = [1 2 3];
params.excit_coefs = [0 5 4 3 1 -1 -2 -2 -1 0 0 0 0 0 0 0 0 0 0 0;
                      1 1 1 1 1  1  1  1  1 1 1 1 1 1 1 1 1 1 1 1];
params.inhib_coefs = [1 1 1 1 1  1  1  1  1 1 1 1 1 1 1 1 1 1 1 1;
                      0 4 3 2 0 -1 -1 -1 -0 0 0 0 0 0 0 0 0 0 0 0];                  
params.filt_dim_num = 2;
params.pretty_name = 'FIR Filters (Pair, Inhib/Excit)';
params.editable_fields = {'num_coefs', 'num_filts', 'excit_sigmoid', ...
                          'inhib_sigmoid'};

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
if ~isequal([params.num_filts params.num_coefs], size(params.inhib_coefs))
    params.excit_coefs = zeros(params.num_filts, params_num_coefs);
    params.inhib_coefs = zeros(params.num_filts, params_num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function pred = do_fir_filtering(parm, x)
    % Apply the FIR filters to corresponding downsampled stimuli to get the model prediction
    % Since it is linear, the prediction is just the sum of the filters
    % We assume that there are no second order terms combining elements of both filters
    inhib_preds = [];
    excit_preds = [];
    for filt_idx = 1:parm.num_filts
        xs = squeeze(x(:,:,filt_idx));
        excit_preds = cat(3, excit_preds, filter(parm.excit_coefs(filt_idx,:), 1, ...
                                       xs, [], parm.filt_dim_num));
        inhib_preds = cat(3, inhib_preds, filter(parm.inhib_coefs(filt_idx,:), 1, ...
                                       xs, [], parm.filt_dim_num));
    end
    % Ensure that the excitory contribution stays positive.
    excit_pred = squeeze(sum(excit_preds, 3));
    excit_pred(le(excit_pred, 0)) = 0;
    
    % Ensure that the inhibitory contribution stays positive
    inhib_pred = squeeze(sum(inhib_preds, 3));    
    inhib_pred(le(inhib_pred, 0)) = 0;
    
    % The predicted response is one minus the other.
    pred =  excit_pred - inhib_pred;
    pred(le(pred, 0)) = 0; % The prediction ALSO cannot go <0
end

function plot_fir_coefs(parm)
    for filt_idx = 1:parm.num_filts
        stem([1:parm.num_coefs], parm.excit_coefs(filt_idx,:), pickcolor((filt_idx-1)*2-1));
        stem([1:parm.num_coefs], parm.inhib_coefs(filt_idx,:), pickcolor((filt_idx-1)*2));
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