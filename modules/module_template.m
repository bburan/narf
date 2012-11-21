function m = module_template(args)
% A template module creation function.
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @module_template;
m.name = 'module_template';
m.fn = @do_elliptic_filter;
m.pretty_name = 'Template Module';
m.editable_fields = {'my_field'};
m.plot_fns = {'Special Plot', @do_plot_special};
m.isready_pred = @module_isready;

% Module fields that are specific to THIS MODULE
m.my_field = 3;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Values computed below this point are not directly editable. Instead, 
% they are computed from the above values in a deterministic manner.

% Finally, define the 'methods' of this module, as if it were a class
function x = do_function(stack, x)
    mdl = stack{end};   % Our module is always the 'top' of the stack

    % Do your modification x here
    x = x + 1;
    
end

function do_plot_special(mdl)    
    plot(1:0.1:10, sin(1:0.1:10));
end

function isready = module_isready(mdl, x)
    isready = false; 
end

end