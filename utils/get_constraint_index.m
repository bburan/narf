function [idx, con] = get_constraint_index(cell_cons, var_name)
% get_constraint_index(cons)
%
% Given constraints \cell_cons\ and a variable name \var\ this function
% returns the cell index and cell content of the matching constraint
%
% ARGUMENTS:
%    cell_cons       a cell array of constraints struct
%    var_name        a string representing the variable name
%
% RETURNS:
%    idx             the cell index
%    con             the constraint struct

idx = 0;
con = struct();

if ~isfield(cell_cons{1},'var')
    error('get_constraint_index must be called with a constraint array!\n');
end
for ii = 1:numel(cell_cons)
    if strcmp(cell_cons{ii}.var, var_name)
        idx = ii;
        con = cell_cons{ii};
        return
    end
end

fprintf('Could not find the constraint...\n');