function generic_data_table(mytable, mystruct, myfields)
% Since data tables are updated in pretty much the same way everywhere, I
% decided to abstract the updating process to avoid code repetition.
l = length(myfields);
c = cell(l,2);
for i = 1:l
    if ~isfield(mystruct, myfields{i})
        log_err('Could not find field: %s', myfields{i});
    end
    c{i,1} = myfields{i};
    c{i,2} = repl_write(mystruct.(myfields{i})); % Ensure data becomes a str
end
set(mytable, 'Data', c);
drawnow;
