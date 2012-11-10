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
%                  params.n_filts    (Number of filters)
%                  params.coefs      (Filter coefficients cell array)
%                  params.freq_resp_plot_fn (Function to plot frq response)
%  
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

% TODO: IFF not already defined....
% Compute the filter coefficients in a cell array 
n_filts = length(params.low_freqs);
params.coefs={}; 
for i = 1:n_filts
    [B,A] = ellip(params.order, 0.5, params.stop_dB, ...
       [(params.low_freqs(i)/params.sampfs)*2,...
        (params.high_freqs(i)/params.sampfs)*2]);   
    params.coefs{i} = {B,A};
end 

% Close over the graphing function below. 
params.freq_resp_plot_fn = @do_plot_elliptic_bandpass_filter_bank_frq_resp;
        
% Return a newly defined function
fn = @do_elliptic_filter;

% This inner function closes over the params variable.
function filtered_x = do_elliptic_filter(x)
% If we didn't get an x, return the param struct of this filter
if (nargin < 1) filtered_x = params; return; end

% TODO: Check the structure of X to ensure its validity





% LEFT OFF HERE: For some reason PARAMS is not showing up here!




filtered_x=[];
for i = 1:n_filts
    tmp = filter(params.coefs{i}{1}, params.coefs{i}{2}, x,[],2);
    filtered_x = cat(3, filtered_x, abs(tmp));
end

% This inner function closes over the 
function do_plot_elliptic_bandpass_filter_bank_frq_resp(ah)
% Accepts an axes handle to plot on.
for filt_idx = 1:n_filts
    ww = 0:(pi/1000):pi;
    H = freqz(params.coefs{filt_idx}{1}, params.coefs{filt_idx}{2}, ww);
    loglog(ww, abs(H), pickcolor(filt_idx));
    setAxisLabelCallback(gca, @(f) (f*SAMPFREQ/(3.14*2)), 'X');
    axis tight;
end

