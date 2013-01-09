function [stack, xxx] = load_model_stack(filename)
    global STACK XXX;
    vars = load(filename, 'stack', 'xxx');
    stack = vars.stack;
    XXX = {};
    XXX{1} = vars.xxx;
    
    % TODO: This is probably not needed but I am superstitious
    for ii = 1:length(stack)            
        if isfield(stack{ii}, 'gh')
            stack{ii} = rmfield(stack{ii}, 'gh');
        end
        if isfield(stack{ii}, 'plot_gui')
            stack{ii} = rmfield(stack{ii}, 'plot_gui');
        end
            
    end
    
end