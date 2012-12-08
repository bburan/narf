function module = generic_model_data_table_update(hDataTable, module)
% Returns an updated copy of a module using values found in hDataTable

s = extract_field_val_pairs(hDataTable, 2, 3);

% If there was an error extracting the fields, reset the GUI
if isempty(s)
   generic_checkbox_data_table(hDataTable, module, module.editable_fields);
   return;
end

for fs = fieldnames(s)', fs=fs{1};
    module.(fs) = s.(fs);
end

module.fit_fields = extract_checked_fields(hDataTable, 1, 2);

module = module.mdl(module);  % Give the module a chance to run its own code


 