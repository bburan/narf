function v = flatten_field(dat, sfs, field)
% Flattens matrices under 'dat.(sf).(field)' into a single long row vector. 
%
% Uses the default order of fields according to matlab (first indexes by
% the row, then by the column, then on to the higher dimensions).
%
% Works on any size matrix

% Count number of elements for each 'field' entry
lens = zeros(length(sfs), 1);
for ii = 1:length(sfs),
    lens(ii) = numel(dat.(sfs{ii}).(field));
end

% Create and fill the row vector
v = zeros(sum(lens), 1);
jj = 1;
for ii = 1:length(sfs),
    v(jj:jj+lens(ii)-1) = reshape(dat.(sfs{ii}).(field), [], 1);
end