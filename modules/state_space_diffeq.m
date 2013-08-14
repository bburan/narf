function m = state_space_diffeq(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @state_space_diffeq;
m.name = 'state_space_diffeq';
m.fn = @do_state_space_diffeq;
m.pretty_name = 'State Space Diff Eq.';
m.editable_fields = {'order', 'element_fn', 'coefs',  ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.n_states = 2; 
m.A_order  = 2;
m.A_fn     = @polyval;
m.B_order  = 1;
m.B_fn     = @times;
m.C_order  = 1;
m.C_fn     = @polyval;
m.D_order  = 1;
m.D_fn     = 1;

m.A_coefs = zeros(m.num_coefs, m.num_dims);
m.B_coefs = zeros(m.num_coefs, m.num_dims);
m.C_coefs = zeros(m.num_coefs, m.num_dims);
m.D_coefs = zeros(m.num_coefs, m.num_dims);
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';
m.init_fit_sig = 'respavg'; % For initializing coefficients only

% Optional fields
m.is_splittable = true;
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_signal(stack, xxx, m.time, m.output);
m.plot_fns{1}.pretty_name = 'Response vs Time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if ~isequal([m.num_dims m.num_coefs], size(m.coefs))
    m.coefs = zeros(m.num_dims, m.num_coefs);
end


% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_state_space_diffeq(stack, xxx)
    mdl = stack{end};
    x = xxx{end};   % Unfortunately, this doesn't seem to copy deeply
    
    % Apply the FIR filter across every stimfile
    for sf = fieldnames(x.dat)', sf=sf{1};
        
        % Compute the size of the filter
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         
         if ~isequal(C, mdl.num_dims)
            error('Dimensions of (mdl.input) don''t match channel count.');
         end

         tmp = zeros(T, S, C);        
         for s = 1:S
             for c = 1:C,
                 % Find proper initial conditions for the filter
                 [~, Zf] = filter(mdl.coefs(c,:)', [1], ...
                                  ones(length(mdl.coefs(c,:)') * 2, 1) .* x.dat.(sf).(mdl.input)(1, s, c));
                              
                 tmp(:, s, c) = filter(mdl.coefs(c,:)', [1], ...
                                       x.dat.(sf).(mdl.input)(:, s, c), ...
                                       Zf);
             end
         end
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = squeeze(sum(tmp, 3)) + mdl.baseline; 
    end
end
end