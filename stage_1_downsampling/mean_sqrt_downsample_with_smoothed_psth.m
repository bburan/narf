function ret = mean_sqrt_downsample_with_smoothed_psth(args)
% Uses the mean of the sqrts of the samples to downsample the stimulus. 
% Uses a gaussian kernel to 'smooth out' the response during downsampling

% Default Parameters
params = [];
params.raw_freq = 100000;        % Original sampling frequency
params.ds_freq = 1000;           % After downsampling, the new sample freq
params.filt_width_in_sigmas = 3; % Plus/minus. Determines accuracy.
params.sigma_in_secs = 0.002;    % How many sec to 'smooth out' the spikes 
params.ds_dim_num = 2;           % The dimension on which to downsample
params.pretty_name = 'mean(sqrt(abs(stim))), filtfilt(gaussian, resp)';
params.editable_fields = {'raw_freq', 'ds_freq', 'filt_width_in_sigmas', ...
                          'sigma_in_secs'};

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

function [ds_stim, ds_resp] = do_mean_sqrt_smooth_downsample(parm, stim, resp)
    % Make a gaussian kernel
    w = 2 * parm.sigma_in_secs * parm.filt_width_in_sigmas * parm.raw_freq; 
    kern = fspecial('gaussian',[1 w], parm.sigma_in_secs * parm.raw_freq); 
    ratio = ceil(parm.raw_freq / parm.ds_freq);
    ds_stim = conv_fn(sqrt(abs(stim)), parm.ds_dim_num, @mean, ratio, 0);
	
    % Filter resp without affecting phase, and preserving magnitude
    tmp = filter(kern, 1, resp, [], 2);
    tmp2 = sqrt(filter(kern, 1, flipdim(tmp, 2), [], 2));
    
    % Downsample this doubly smoothed filter by just picking nth points
    ds_resp = conv_fn(tmp2, parm.ds_dim_num, @(x) (x(1)), ratio, 0);
    
    disp('dbg');
end

% ------------------------------------------------------------------------
% Put the instance methods in params struct
params.downsamp_fn = @do_mean_sqrt_smooth_downsample;

ret = params;
end