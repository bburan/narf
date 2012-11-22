function s = extract_checked_fields(mytable, checkbox_col, fieldname_col)
% Return a cell array of fields with checked boxes next to them.
d = get(mytable, 'Data');
[r, c] = size(d);
j = 1;
s = {};
for i = 1:r
    if d{i,checkbox_col}
        s{j} = d{i,fieldname_col};
        j = j+1;
    end
end