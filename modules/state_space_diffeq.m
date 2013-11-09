function m = state_space_diffeq(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @state_space_diffeq;
m.name = 'state_space_diffeq';
m.fn = @do_state_space_diffeq;
m.pretty_name = 'State Space Diff Eq.';
m.editable_fields = {'A', 'B', 'C', 'D', 'delay_B', 'delay_B_amount', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.A = [-0.1 -.1; .1 -0.1];
m.B = [0 0; 0 0];
m.C = [1 0];
m.D = [0 0];
m.delay_B        = [0 0; 1 -1];
m.delay_B_amount = 0.5;
m.x_0     = [0 0];
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_single_default_output;
m.plot_fns{1}.pretty_name = 'Output vs Time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_state_space_diffeq(mdl, x, stack, xxx)
    
    % Create a delay term matrix just to delay the inputs slightly.
    delayterms = struct('delay', abs(mdl.delay_B_amount), ...
        'a',[],'b', mdl.delay_B,'c', [],'d', []); 

    try
        sys=delayss(mdl.A, mdl.B, mdl.C,mdl.D, delayterms);
    catch
        keyboard
    end
    
    % Apply the FIR filter across every stimfile
    for sf = fieldnames(x.dat)', sf=sf{1};
        
        % Compute the size of the filter
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         
         if ~isequal(C, size(mdl.A, 1))
            error('Dimensions of A don''t match channel count.');
         end
             
         tmp = zeros(T, S, 1);        
         
         for s = 1:S
             
             % TODO: FIND x_0 for each model?
             
             tmp(:,s) = lsim(sys, ...
                                squeeze(x.dat.(sf).(mdl.input)(:, s, :)), ...
                                x.dat.(sf).(mdl.time)(:,1), ...
                                mdl.x_0);
             
         end
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = tmp;
    end
end
end