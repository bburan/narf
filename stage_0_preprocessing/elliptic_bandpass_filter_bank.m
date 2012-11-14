function ret = elliptic_bandpass_filter_bank(args)
% This function is a way of creating pseudo-object-oriented blobs of
% functionality. Basically, a zero-argument version of this code returns
% the default elliptic bandpass filter object. If you pass it a single
% argument, the 'params' struct, a new object will be created using those
% values and filling in the defaults where necessary. 
%
% The reason this is pseudo-object oriented is that you can give the return
% 'params' struct function handles which will act like instance methods,
% without all the hassle of actually learning OOP on matlab. 
% 
% This is more like programming to an interface, since we don't care what
% exact 'type' of blob it is, just what instance methods it possesses. 
%
% Methods that NARF_GUI will use to interact with a blob follow:
%    preproc_fn(x)              where x is a matrix.... (TODO)
%    downsamp_fn(x)             where x is ... (TODO)
%    model_fn(x)                where x is ... (TODO)
%    
% Other interface functions include:
%    freq_resp_plot_fn(params)  
%    
% ----------------------------------------------------------------------
% Parameters specific to ELLIPTIC_BANDPASS_FILTER_BANK
% low_freqs  (Bottoms of the bandpass filters)
% high_freqs (Tops of the bandpass filters)
% order      (Order of the filter)
% stop_dB    (Minimum stopband attenuation)
% sampfs     (Sampling frequency in Hz)
%
% Also, few parameters are not modifiable and will always be automanaged:
% coefs            (Filter coefficients cell array)
% pretty_name      (User readable string)
% editable_fields  (A cell array of editable parameter field names)

% Default Parameters
params = [];
params.low_freqs = [2000 20000];
params.high_freqs = [4000 27000];    
params.order = 4; 
params.sampfs = 100000; 
params.stop_dB = 50; 
params.pretty_name = 'Elliptic Bandpass Filter Bank';
params.editable_fields = {'low_freqs', 'high_freqs', 'order', 'sampfs', 'stop_dB'};

% Overwrite the defaults with the arguments iff they are in editable_fields
if nargin == 1
    fns = fieldnames(args);
    for idx = 1:length(fns);
        if any(ismember(params.editable_fields, fns{idx}))
            params.(fns{idx}) = args.(fns{idx});
        end
    end
end

% Below values are computed using above values and are not directly settable
params.coefs={}; 
for i = 1:length(params.low_freqs);
    [B,A] = ellip(params.order, 0.5, params.stop_dB, ...
                  [(params.low_freqs(i)/params.sampfs)*2,...
                   (params.high_freqs(i)/params.sampfs)*2]);   
    params.coefs{i} = {B,A};
end 

% ------------------------------------------------------------------------
% INSTANCE METHODS FOLLOW

function filtered_x = do_elliptic_filter(parm, x)
    % If we didn't get an x, return the param struct of this filter
    if (nargin ~= 2) error('WRONG NUMBER OF ARGS STUPID!'); end
    % TODO: Check the structure of x to ensure its validity
    filtered_x=[];
    for i = 1:length(parm.low_freqs)
        tmp = filter(parm.coefs{i}{1}, parm.coefs{i}{2}, x,[],2);
        filtered_x = cat(3, filtered_x, tmp); 
    end
end

function do_plot_elliptic_bandpass_filter_bank_frq_resp(parm)    
    if (nargin ~= 1) error('WRONG NUMBER OF ARGS STUPID!'); end
    for filt_idx = 1:length(parm.low_freqs)
        ww = 0:(pi/1000):pi;
        H = freqz(parm.coefs{filt_idx}{1}, parm.coefs{filt_idx}{2}, ww);
        loglog(ww, abs(H), pickcolor(filt_idx));
        setAxisLabelCallback(gca, @(f) (f*parm.sampfs/(3.14*2)), 'X');
        axis tight;
    end 
end

% TODO: Make this a validation function which can be called before the
% actual filtering process. 
% Check that the merged parameters are valid
% if ~(all(isfield(params, {'high_freqs', 'low_freqs', 'order', ...
%         'stop_dB', 'sampfs'})) & ...
%      isvector(params.low_freqs) & isvector(params.high_freqs) & ...
%      length(params.low_freqs) == length(params.high_freqs))
%     error('Elliptic bandpass frequency settings are invalid or missing.');
% end

% ------------------------------------------------------------------------
% Put the instance methods in the params struct
params.freq_resp_plot_fn = @do_plot_elliptic_bandpass_filter_bank_frq_resp;
params.preproc_fn = @do_elliptic_filter;

ret = params;

end
