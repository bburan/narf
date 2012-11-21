function mods = scan_directory_for_modules(scanpath)
% Scan a directory for executable module objects
% If any exception occurs when loading a module, it will be dropped.

fprintf('Scanning dir for modules: %s\n', scanpath);
files = dir(fullfile(scanpath, '*.m'));
mods = [];
for i = 1:length(files)
    f = files(i).name;
    try
        name = f(1:end-2);    % File name minus the '.m' at end
        mod = str2func(name);  % Make an executable function
        m = mod();             % Zero args means "ask for defaults"
        % Check that all required interface fields are defined
        if ~(all(isfield(m, ...
                {'mdl', ...               % Module function handle
                 'name', ...              % Module function name as string
                 'pretty_name',  ...      % Module name displayed to users
                 'fn', ...                % Function handle to execute
                 'editable_fields', ...   % Which params are editable
                 'plot_fns', ...          % Plot_fn param structs
                 'isready_pred'})) & ... % Fn's executability predicate
                  isequal(name, m.name))
            error('Modules must define all of the following fields: mdl, fn, name, pretty_name, editable_fields, plot_mods, isready_pred');
        % TODO: Check the structure of plot_fns as well!
        else     
           fprintf('Found ''%s''\n', f);
           mods.(name) = m;  % Index module under its module name
        end
    catch ME
        fprintf('Skipping %s due to error: %s\n', f, ME.message);
    end
end
