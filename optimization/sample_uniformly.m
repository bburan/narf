function samples = sample_uniformly(a, b, n)
% Samples n times between lower limit vector a, upper limit vector b,
% with a uniform probability density function between the two. 
% Returns a matrix where each row is a useful sample vector.
%
% Oct 12, 2012. Ivar Thorson.

if size(a) ~= size(b)
    error('Sizes of the minimum vector a and maximum vector b must match');
end

if ~(all(b > a))
    error('Min vector a must be strictly smaller than max vector b.');
end

[nr nc] = size(a);

% To make this work regardless of whether a,b were row or col vectors
if (nr < nc) 
    samples = repmat(a, n, 1) + repmat((b-a), n, 1) .* rand(n, length(a));
else
    samples = repmat(a', n, 1) + repmat((b-a)', n, 1) .* rand(n, length(a));
end