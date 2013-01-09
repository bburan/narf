function load_model_stack(filename)
    global STACK XXX;
    vars = load(filename, 'stack', 'xxx');
    
    XXX = {};
    XXX{1} = vars.xxx;
    stack = vars.stack;
    
    % Connect the old GUI to the new model, if the old gui exists
    for ii = 1:length(stack)
        if isfield(stack{ii}, 'gh')
            STACK{ii}.gh = stack{ii}.gh;
        end
        if isfield(stack{ii}, 'plot_gui')
            STACK{ii}.plot_gui = stack{ii}.plot_gui;
        end
    end
    
    % Swap in the new stack
    STACK = vars.stack;
end