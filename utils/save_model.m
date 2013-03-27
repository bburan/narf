function save_model(filename, stack, xxx, meta)
% Save the global STACK, XXX, and META data structures to disk

if nargin ~= 4
    error('save_model() needs exactly 4 arguments');
end 

% Strip off any GUI handles and plot GUI handles
for ii = 1:length(stack)
    if isfield(stack{ii}, 'gh')
        stack{ii} = rmfield(stack{ii}, 'gh');
    end
    if isfield(stack{ii}, 'plot_gui')
        stack{ii} = rmfield(stack{ii}, 'plot_gui');
    end
end

% Save to disk
xxx = xxx{1};
save(filename, 'stack', 'xxx', 'meta');

% Change the file to be read-only so you don't accidentally alter it later
unix(['chmod 666 ' filename]);