function v = flatten_field(dat, sfs, field)
% v = flatten_field(dat, sfs, field)
%
% Concatenates each matricies found under 'dat.(sf).(field)', where sf is
% each of the elements of cell array sfs. Returns a single long vector.
% Intended to be used for quickly creating a 'resp' or 'stim' type vector.
%
% ARGUMENTS:
%    dat    A structure
%    sfs    Structure fieldnames to be used
%    field  The field to be extracted
%
% RETURNS:
%    v      A single long vector containing all the values
%
% NOTE: Relies on the default order of fields in each matrix being
% organized such that the first dimenson is time. 

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
    jj=jj+lens(ii);
end