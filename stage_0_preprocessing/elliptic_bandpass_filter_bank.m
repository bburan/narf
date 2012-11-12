function fn = elliptic_bandpass_filter_bank(params)
% Return a function  FN that applies a sequence of elliptic bandpass
% filters to its input. FN accepts two arguments:
%    1. A 1xN vector of the stimuli vector
%    2. The sampling frequency 
% 
% FN will also return a single matrix, where the rows are the time stimulus
% filtered by the F'th filter. The number of filters is determined by the
% params vector given in this function.
%
% INPUT: 
%    params    A struct with at lesat these fields:
%                  params.low_freqs  (Bottoms of the bandpass filters)
%                  params.high_freqs (Tops of the bandpass filters)
%                  params.order      (Order of the filter)
%                  params.stop_dB    (Minimum stopband attenuation)
%                  params.sampfs     (Sampling frequency in Hz)
%              Optionally, it may also have these fields:
%                  params.coefs      (Filter coefficients cell array)
%                  params.freq_resp_plot_fn (Function to plot frq response)
%  

% TODO: Make this initialization work even if we only get SOME params but
% not all of them. 

% If no params are given get the defaults.
if (nargin < 1)
    params = [];
    % TODO: Check that there is info in the CFD files about SPN band sizes,
    % and initialize freq bands accordingly
    params.low_freqs = [2000 20000];
    params.high_freqs = [4000 27000];
    params.order = 4;
    params.stop_dB = 50;
    params.sampfs = 100000;
    params.pretty_name = 'Elliptic Bandpass Filter Bank';
    params.freq_resp_plot_fn = @do_plot_elliptic_bandpass_filter_bank_frq_resp;
    params.editable_fields = {'low_freqs', 'high_freqs', 'order', ...
                              'stop_dB'};
    fn = params;
    return;
end

% Check that the parameters are valid
if ~(all(isfield(params, {'high_freqs', 'low_freqs', 'order', ...
        'stop_dB', 'sampfs'})) & ...
     isvector(params.low_freqs) & isvector(params.high_freqs) & ...
     length(params.low_freqs) == length(params.high_freqs))
    error('Elliptic bandpass frequency settings are invalid or missing.');
end

% Compute the filter coefficients in a cell array 
n_filts = length(params.low_freqs);
params.coefs={}; 
for i = 1:n_filts
    [B,A] = ellip(params.order, 0.5, params.stop_dB, ...
       [(params.low_freqs(i)/params.sampfs)*2,...
        (params.high_freqs(i)/params.sampfs)*2]);   
    params.coefs{i} = {B,A};
end 

% ------------------------------------------------------------------------
% DEFINE TWO INNER FUNCTIONS...who needs objects when you can use closures!

% BEGIN INNER FUNCTION: 
function filtered_x = do_elliptic_filter(x)
    % If we didn't get an x, return the param struct of this filter
    if (nargin < 1) filtered_x = params; return; end
    % TODO: Check the structure of x to ensure its validity
    filtered_x=[];
    for i = 1:n_filts
        tmp = filter(params.coefs{i}{1}, params.coefs{i}{2}, x,[],2);
        filtered_x = cat(3, filtered_x, tmp); % TODO: Do abs here or not?
    end
end % END INNER FUNCTION

% BEGIN INNER FUNCTION (plot the response)
function do_plot_elliptic_bandpass_filter_bank_frq_resp(ah, params)
    % AH: Axes handle to plot on.
    % PARAMS: Should be the same as above.
    for filt_idx = 1:n_filts
        ww = 0:(pi/1000):pi;
        H = freqz(params.coefs{filt_idx}{1}, params.coefs{filt_idx}{2}, ww);
        loglog(ww, abs(H), pickcolor(filt_idx));
        setAxisLabelCallback(gca, @(f) (f*params.sampfs/(3.14*2)), 'X');
        axis tight;
    end 
end % END THE INNER FUNCTION
% ------------------------------------------------------------------------

% Return newly defined inner function
fn = @(arg) do_elliptic_filter(arg);

end



