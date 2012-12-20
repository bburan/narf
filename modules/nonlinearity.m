function m = nonlinearity(args)
% Applies a general nonlinear function to the input. 
%
% The nonlinear function (NLFN) must be a function of two arguments:
%
%    NLFN(PHI, X) -> Y
%
% where
%    PHI is NLFN's parameter vector, such as polynomial coefficients
%    X is the matrix that the function works on
%    Y is the output of NLFN
%    
% Good choices for NLFN are 
%    polyval         (See matlab documentation)
%    sigmoid         (See utils/ directory)
%    exponential     (See utils/ directory)
%
% You may of course define your own functions. 

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonlinearity;
m.name = 'nonlinearity';
m.fn = @do_nonlinearity;
m.pretty_name = 'Generic Nonlinearity';
m.editable_fields = {'input', 'time', 'output', 'nlfn', 'phi'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time = 'stim_time';
m.output = 'stim';
m.nlfn = @polyval;
m.phi = [1 0];   % Default is to pass the signal through untouched

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_output_vs_time(stack, xxx, m.time, m.output);
m.plot_fns{1}.pretty_name = 'Output vs Time';

m.plot_fns{2}.fn = @(stack, xxx) do_plot_nonlinearity(stack, xxx, stack{end}.input, @(x) stack{end}.nlfn(stack{end}.phi, x), false);
m.plot_fns{2}.pretty_name = 'Nonlinearity';

m.plot_fns{3}.fn = @(stack, xxx) do_plot_nonlinearity(stack, xxx, stack{end}.input, @(x) stack{end}.nlfn(stack{end}.phi, x), true);
m.plot_fns{3}.pretty_name = 'Nonlinearity + Histogram';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        y = zeros(T, S, C);
       
        % TODO: If a scalar-valued function, use this
        %y = arrayfun(@(in) mdl.nlfn(mdl.phi, in), x.dat.(sf).(mdl.input));
        
        % Otherwise use the much faster vector valued functions
        y = mdl.nlfn(mdl.phi, x.dat.(sf).(mdl.input));
        
        x.dat.(sf).(mdl.output) = y;
    end
end

end