function m = reparameterize(args)
% Allows you to generate a new parameterization for a different module

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @reparameterize;
m.name = 'reparameterize';
m.fn = @do_reparameterize;
m.pretty_name = 'Reparameterize';
m.editable_fields = {'params', 'param_fn', 'output_field'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.params       = [pi/6 3];
m.param_fn     = @(phi) sin(phi(1)) * phi(2);
m.output_field = 'output_field';

% Optional fields
m.plot_fns = {};

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_reparameterize(mdl, x, stack, xxx)    
    % Super duper simple: 
    x.(mdl.output_field) = mdl.param_fn(mdl.params);
end

end