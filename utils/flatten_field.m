function v = flatten_field(dat, field)
% Flattens matrices under 'dat.*.(field)' into a single long row vector. 
%
% Uses the default order of fields according to matlab (first indexes by
% the row, then by the column, then on to the higher dimensions).
%
% Works on any size matrix

% Count number of elements for each 'field' entry
fns = fieldnames(dat);
lens = zeros(length(fns));
for ii = length(fns),
    lens(ii) = numel(dat.(fns{ii}).(field));
end

% Create and fill the row vector
v = zeros(1, sum(lens));
jj = 1;
for ii = length(fns),
    v(jj:jj+lens{ii}) = reshape(dat.(fns{ii}).(field), [], 1);
end