function delete_module_gui(idx)
% delete_module_gui()
%
% Deletes any gui handles, widgets, etc from the module at STACK index
% number idx. 
%
% No arguments or return values. Use purely for side effects. 

global STACK;

if isfield(STACK{idx}, 'gh')
    try
        delete(STACK{idx}.gh.plot_axes);
        delete(STACK{idx}.gh.plot_popup);
        delete(STACK{idx}.gh.plot_panel);
        delete(STACK{idx}.gh.fn_apply);
        delete(STACK{idx}.gh.fn_table);
        delete(STACK{idx}.gh.fn_popup);
        delete(STACK{idx}.gh.fn_panel);
    catch
        % Do nothing if a delete failed
    end
    STACK{idx} = rmfield(STACK{idx}, 'gh');
end
if isfield(STACK{idx}, 'plot_gui')
    STACK{idx} = rmfield(STACK{idx}, 'plot_gui');
end
