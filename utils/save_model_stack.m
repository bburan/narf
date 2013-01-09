function save_model_stack(filename, stack, xxx)
    % Saves the model stack and XXX{1} initialization info
   
    for ii = 1:length(stack)
        % Strip off any GUI handles
        if isfield(stack{ii}, 'gh')
            stack{ii} = rmfield(stack{ii}, 'gh');
        end
        % Strip off any plot GUI handles
        if isfield(stack{ii}, 'plot_gui')
            stack{ii} = rmfield(stack{ii}, 'plot_gui');
        end
    end
    
    % Save the stack and initialization
    xxx = xxx{1};
    save(filename, 'stack', 'xxx');
    
end