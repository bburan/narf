function generic_checkbox_data_table(mytable, mystruct, myfields)
% Since data tables are updated in pretty much the same way everywhere, I
% decided to abstract the updating process to avoid code repetition.
l = length(myfields);
c = cell(l,3);
for i = 1:l
    if ~isfield(mystruct, myfields{i})
        error('Could not find field: %s', myfields{i});
    end
    if isfield(mystruct, 'fit_fields')
        c{i,1} = isfield(mystruct.fit_fields, myfields{i});
    else
        c{i,1} = false;
    end
    c{i,2} = myfields{i};
    c{i,3} = repl_write(mystruct.(myfields{i})); % Ensure data becomes a str
end
set(mytable, 'ColumnName', {'Fit?', 'Field', 'Value'});
set(mytable, 'ColumnEditable', [true false true]);
set(mytable, 'ColumnWidth', {25 150 100});
set(mytable, 'RowName', {});
set(mytable, 'Data', c);

drawnow;
