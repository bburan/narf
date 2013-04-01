function z = iconv(x, f)
% z = iconv(x, f)
%
% Ivar's boundary-preserving convolution function. Effectively similar to
% using CONV(X, F, 'same') except that it preserves the edges instead of
% smoothing them towards zero. 
%
% ARGUMENTS:
%    x    The signal vector you want to convolve
%    f    A filter to convolve with
%
% RETURNS:
%    z    The convolved signal. It will have the same length as f. 
%

z = zeros(size(x));
lx = length(x);
lf = length(f);
tmp = conv(x, f, 'full');
lt = length(tmp);

% Copy the center part of the filter
for ii = 1:lx
    z(ii) = tmp(ii+lf/2);
end

% "Double up" the ends of the filter
for ii = 1:lf/2
    z(ii) = (z(ii) + tmp((lf/2) - ii + 1)); % Left side
    z(end-ii+1) = (z(end-ii+1) + tmp(end - (lf/2) + ii)); % Right side
end

