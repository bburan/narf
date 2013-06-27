function ret = sort_by_field(cellarray, field)
% ret = sort_by_field(cellarray, field)
% 
% Returns a sorted copy of CELLARRAY with structs sorted by the value 
% stored under fieldname FIELD.
%
% ARGUMENTS:
%    cellarray   A cell array of structs.
%    field       A string indicating the fieldname to sort by. 
%
% RETURNS:
%    ret         A sorted copy of cellarray.

ff = getfieldforcellarray(cellarray, field);
[~, idxs] = sort(cell2mat(ff));
ret = cellarray(idxs);