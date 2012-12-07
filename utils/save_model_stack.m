function save_model_stack(filename, stack)
   
    for ii = 1:length(stack)
        % Strip off any GUI handles
        if isfield(stack{ii}, 'gh')
            rmfield(stack{ii}, 'gh');
        end
        % Strip off any plot GUI handles
        if isfield(stack{ii}, 'plot_gui')
            rmfield(stack{ii}, 'plot_gui');
        end
    end
    
    % Save the stack
    save(filename, 'stack');
    
end