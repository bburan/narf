function s = extract_field_val_pairs(mytable, fieldname_col, value_col)
% Return a new struct extracted from two columns
d = get(mytable, 'Data');
[r, c] = size(d);
if (c < fieldname_col | c < value_col | fieldname_col < 1 |  value_col < 1)
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