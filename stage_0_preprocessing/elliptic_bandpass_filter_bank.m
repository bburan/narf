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
%    params    A struct with all of these fields:
%                  params.low_freqs  (Bottoms of the bandpass filters)
%                  params.high_freqs (Tops of the bandpass filters)
%                  params.order      (Order of the filter.)
%                  params.stop_dB    (Minimum stopband attenuation.)

% If no params are given get the defaults. 
if (nargin < 1)
    params = [];
    % TODO: Check that there is info in the CFD files about SPN band sizes,
    % and initialize freq bands accordingly
    params.low_freqs = [2000 20000];
    params.high_freqs = [4000 27000];
    params.order = 4;
    params.stop_dB = 50;    
end

% Check that the parameters are valid
if ~(all(isfield(params, ['high_freqs', 'low_freqs', 'order', 'stop_dB'])) & ...
     isvector(params.low_freqs) & isvector(params.high_freqs) & ...
     length(params.low_freqs) == length(params.high_freqs))
    error('Elliptic bandpass frequency settings are invalid or missing.');
    return;
end

% Closed-over variable #1: Number of filters
n_filts = length(params.low_freqs);

% Closed-over variable #2: Filter coefficients in a cell array.
coefs={}; 
for i = 1:n_filts
    [B,A] = ellip(params.order, 0.5, params.stop_dB, ...
       [(params.low_freqs(i)/SAMPFREQ)*2, (params.high_freqs(i)/SAMPFREQ)*2]);   
    coefs{i} = {B,A};
end

% Return a fn which closes over the above.
fn = @do_elliptic_filter;

% The inner function which closes over the above two variables
function filtered_x = do_elliptic_filter(x)
% TODO: Check the structure of X to ensure its validity
filtered_x=[];
for i = 1:n_filts
    tmp = filter(coefs{i}{1}, coefs{i}{2}, x,[],2);
    filtered_x = cat(3, filtered_x, abs(tmp));
end
