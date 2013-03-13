function delete_all_module_guis()
global STACK;
    for ii = 1:length(STACK)
        if isfield(STACK{ii}, 'gh')
            try
                delete(STACK{ii}.gh.plot_axes);
                delete(STACK{ii}.gh.plot_popup);
                delete(STACK{ii}.gh.plot_panel);
                delete(STACK{ii}.gh.fn_apply);
                delete(STACK{ii}.gh.fn_table);
                delete(STACK{ii}.gh.fn_popup);
                delete(STACK{ii}.gh.fn_panel);
            catch
                % Do nothing if a delete failed
            end
            STACK{ii} = rmfield(STACK{ii}, 'gh');
        end
        if isfield(STACK{ii}, 'plot_gui')
            STACK{ii} = rmfield(STACK{ii}, 'plot_gui');
        end
    end
end
