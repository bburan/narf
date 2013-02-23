function z = iconv(x, f)
% iconv(x, f) : Ivar's boundary-preserving convolution function.
% Probably just a hack, but it's a quick way to smooth. 
%
% ARGUMENTS:
%    x, the signal vector you want to convolve with
%    f, a filter to convolve with
%
% The returned signal vector is the same length as x, but should have edges
% that are twice as tall as would have been created with conv(x, f, 'same')

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
    
    % figure;
    % plot(1:lt, tmp, 'k-', lf/2:lt-lf/2, x, 'r-');
    % plot(1:lx, x, 'r-', 1:lx, z);
    
end    
