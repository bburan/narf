function s = extract_field_val_pairs(mytable, fieldname_col, value_col)
% s = extract_field_val_pairs(mytable, fieldname_col, value_col)
% 
% Return a new struct s whose fieldnames are extracted from one column and
% whose values are extracted from another. The table is assumed to hold
% strings which will be processed with matlab EVAL().
%
% ARGUMENTS:
%    mytable         A handle to a UITABLE object
%    fieldname_col   Index of the column containing the field names
%    value_col       Index of the column containing the values.
%                    Values are assumed to be strings which can be
%                    converted into matlab values via EVAL().
%
% RETURNS:
%    s               A struct c
%

d = get(mytable, 'Data');
[r, c] = size(d);
if (c < fieldname_col || c < value_col || fieldname_col < 1 ||  value_col < 1)
    err('Column index number is outside the data table''s range.');
end
s = {};
for i = 1:r
    try 
        s.(d{i,fieldname_col}) = eval(d{i,value_col});
    catch
       % If there was a problem, return an empty struct
       s = {};
       return
    end
end