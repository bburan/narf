function ret = fir_filter_volterra(args)
% A 2D FIR Filter is created for each input filter dimension.
% The total number of filter coefficients is: num_coefs*(num_filts)^2
% (Sorry about the naming, num_coefs is the number of coefficients along
% the time dimension of the filter. num_coefs=10 remembers the last 10 
% samples, for example)

% Default Parameters
params = [];
params.num_coefs = 20;
params.num_filts = 2;
params.coefs = {};
for i = 1:params.num_coefs
    params.coefs{i} = [0 5 4 3 1 -1 -2 -2 -1 0 0 0 0 0 0 0 0 0 0 0;
                       1 1 1 1 1  1  1  1  1 1 1 1 1 1 1 1 1 1 1 1];
end
params.filt_dim_num = 2;
params.pretty_name = 'FIR Filters (2nd Order Volterra)';
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
% if ~isequal([params.num_filts params.num_coefs], size(params.coefs))
%     params.coefs = zeros(params.num_filts, params_num_coefs);
% end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function pred = do_fir_volterra_filter(parm, x)
    % Applies multiple volterra filters to the input vector, then sums the
    % responses and returns the prediction
     %preds = [];
%      for f1_idx = 1:parm.num_filts
%          % 
%          for f2_idx = 1:parm.num_filts
%              p
%          end
%          % Combine those 
%          preds = cat(3, preds, filter(parm.coefs(filt_idx,:), 1, ...
%                                         squeeze(x(:,:,filt_idx)), [], ...
%                                         parm.filt_dim_num));
%      end
%      pred = squeeze(sum(preds, 3));
        pred = [];
end

function plot_fir_coefs(parm)
%     for f1_idx = 1:parm.num_filts
%         for f2_idx = 1:parm.num_filts
%             stem([1:parm.num_coefs], parm.coefs(filt_idx,:), pickcolor(filt_idx));
%         end
%     end
%     
    %setAxisLabelCallback(gca, @(t) (t-1), 'X');
    axis tight;
end

% ------------------------------------------------------------------------
% Put the instance methods in params struct
params.model_fn = @do_fir_volterra_filter;
params.plot_coefs_fn = @plot_fir_coefs;
ret = params;
end