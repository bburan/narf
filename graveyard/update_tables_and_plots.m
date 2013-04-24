function update_tables_and_plots()

global STACK;

for ii = 1:length(STACK);
    generic_checkbox_data_table(STACK{ii}.gh.fn_table, STACK{ii}, STACK{ii}.editable_fields); 
    
    hgfeval(get(STACK{ii}.gh.plot_popup, 'Callback'), ...
            STACK{ii}.gh.plot_popup, []);
end
