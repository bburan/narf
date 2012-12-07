function stack = load_model_stack(filename)
    global STACK XXX;
    XXX = XXX(1);
    stack = load(filename, 'stack');
    stack = stack.stack;
%     
%     for ii = 1:length(stack)
%         % Connect to an existing GUI, if one exists
%         if isfield(STACK{ii}, 'gh')
%             stack{ii}.gh = STACK{ii}.gh;
%         end
%         % Strip off any plot GUI handles
%         if isfield(STACK{ii}, 'plot_gui')
%             stack{ii}.plot_gui = STACK{ii}.plot_gui;
%         end
%     end
    
end