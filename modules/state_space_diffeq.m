function m = state_space_diffeq(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @state_space_diffeq;
m.name = 'state_space_diffeq';
m.fn = @do_state_space_diffeq;
m.pretty_name = 'State Space Diff Eq.';
m.editable_fields = {'A', 'B', 'C', 'D', 'delay_B', 'delay_B_amount', 'y_offset', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.A = [-0.1 -.1; .1 -0.1];
m.B = [0 0; 0 0];
m.C = [1 0];
m.D = [0 0];
m.delay_B        = [0 0; 1 -1];
m.delay_B_amount = 20;  % In ms
m.x_0     = [0 0];
m.y_offset = 0;
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_abcd_impulse_response;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_single_default_output;
m.plot_fns{1}.pretty_name = 'Output vs Time';
m.plot_fns{2}.fn = @do_plot_abcd_impulse_response;
m.plot_fns{2}.pretty_name = 'Impulse Response';
m.plot_fns{3}.fn = @do_plot_abcd_step_response;
m.plot_fns{3}.pretty_name = 'Step Response';
m.plot_fns{4}.fn = @do_plot_abcd_bodemag_plot;
m.plot_fns{4}.pretty_name = 'Bode Mag. Plot';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function sys = makesys(mdl)
    % Create a delay term matrix just to delay the inputs slightly.
    delayterms = struct('delay', abs(mdl.delay_B_amount ./ 1000), ...
        'a',[],'b', mdl.delay_B,'c', [],'d', []); 
    sys = delayss(mdl.A, mdl.B, mdl.C,mdl.D, delayterms);
end

function x = do_state_space_diffeq(mdl, x, stack, xxx)
    
    sys = makesys(mdl);
    
    % Apply the FIR filter across every stimfile
    for sf = fieldnames(x.dat)', sf=sf{1};
        
        % Compute the size of the filter
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         
         if ~isequal(C, size(mdl.A, 1))
            error('Dimensions of A don''t match channel count.');
         end
             
         tmp = zeros(T, S, 1);
         
         for s = 1:S
             x0 = mdl.x_0;    % TODO: FIND x_0 for each stim case. 
             t = x.dat.(sf).(mdl.time)(:,1);
             u = squeeze(x.dat.(sf).(mdl.input)(:, s, :));     
             warning off Control:analysis:LsimStartTime;
             tmp(:,s) = lsim(sys, u, t, x0);             
             warning on Control:analysis:LsimStartTime;
         end
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = tmp + mdl.y_offset;;
    end
end

function do_plot_abcd_impulse_response(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    impulse(sys);   
    do_xlabel('Time [s]');
    do_ylabel('Impulse Response');
end

function do_plot_abcd_step_response(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};    
    sys = makesys(mdls{1});
    step(sys);
    do_xlabel('Time [s]');
    do_ylabel('Step Response');
    
end

function do_plot_abcd_bodemag_plot(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    bodemag(sys);
    do_xlabel('Freq [Hz]');
    do_ylabel('Freq Response');
end

end