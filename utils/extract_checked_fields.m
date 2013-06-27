function c = extract_checked_fields(mytable, checkbox_col, fieldname_col)
% c = extract_checked_fields(mytable, checkbox_col, fieldname_col)
% 
% Return a new cell array c of fieldnames of rows which have checked
% checkboxes in in them.
%
% ARGUMENTS:
%    mytable         A handle to a UITABLE object
%    checkbox_col    Index of the column containing checkbox widgets
%    fieldname_col   Index of the column containing the field names
%
% RETURNS:
%    c         Cell array of the fieldnames of rows with checked checkboxes
%
% Return a cell array of fields with checked boxes next to them.

d = get(mytable, 'Data');
r = size(d, 1);
c = {};
for i = 1:r
    if d{i,checkbox_col}
        c{end+1} = d{i,fieldname_col};
    end
end