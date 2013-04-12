function save_model(filepath, stack, xxx, meta)
% save_model(filepath, stack, xxx, meta)
%
% Saves the provided STACK, XXX, and META data structures to disk so that
% they can be loaded again later with load_model();
%
% ARGUMENTS:
%    filepath   The absolute path to where you want to save the model.
%    stack      The STACK data structure you want to save
%    xxx        The XXX data structure, of which only the 1st element will
%               be saved because the rest can always be recreated later. 
%    meta       Metadata about the model, such as its name, git hash tag,
%               fitter name, etc
%
% RETURNS: Nothing

if nargin ~= 4
    error('save_model() needs exactly 4 arguments');
end

% OBSOLETED now that GUI stuff is in a separate structure NARFGUI
% % Strip off any GUI handles and plot GUI handles
% for ii = 1:length(stack)
%     if isfield(stack{ii}, 'gh')
%         stack{ii} = rmfield(stack{ii}, 'gh');
%     end
%     if isfield(stack{ii}, 'plot_gui')
%         stack{ii} = rmfield(stack{ii}, 'plot_gui');
%     end
% end

% Save to disk
xxx = xxx{1};
if exist(filepath,'file'),
    delete(filepath);
end
save(filepath, 'stack', 'xxx', 'meta');

% Change the file to be read-only so you don't accidentally alter it later
% unix(['chmod 444 ' filepath]);
