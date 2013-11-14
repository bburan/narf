function m = pole_zeros(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @pole_zeros;
m.name = 'pole_zeros';
m.fn = @do_pole_zeros;
m.pretty_name = 'Pole/Zeros';
m.editable_fields = {'poles', 'B', 'delay', 'y_offset', 'delay_per_chan', ...
                     'input', 'time', 'output'};  
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.order = 3;
m.poles = [];
m.A = [];
m.B = [];
m.C = [];
m.D = [];
m.delay = [20]; % ms
m.delay_per_chan = false;
m.x_0     = [];
m.y_offset = 0;
m.num_inputs = 0;
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_pz_impulse_response;
m.auto_init = @auto_init_pz;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_single_default_output;
m.plot_fns{1}.pretty_name = 'Output vs Time';
m.plot_fns{2}.fn = @do_plot_pz_impulse_response;
m.plot_fns{2}.pretty_name = 'Impulse Response';
m.plot_fns{3}.fn = @do_plot_pz_step_response;
m.plot_fns{3}.pretty_name = 'Step Response';
m.plot_fns{4}.fn = @do_plot_pz_bodemag_plot;
m.plot_fns{4}.pretty_name = 'Bode Mag. Plot';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function mdl = auto_init_pz(stack, xxx)
    % NOTE: Unlike most plot functions, auto_init functions get a 
    % STACK and XXX which do not yet have this module or module's data
    % added to them.    
    if ~isfield(m, 'fit_fields') 
        return
    end
    
    % Init the coefs to have the right dimensionality
    mdl = m;
    x = xxx{end};
    fns = fieldnames(x.dat);
    sf = fns{1};
    [T, S, C] = size(x.dat.(sf).(mdl.input));               
    
    mdl.num_inputs = C;         
    mdl.poles = ones(mdl.order,1);
    mdl.A = zeros(mdl.order, mdl.order); 
    mdl.B = zeros(mdl.order, mdl.num_inputs + 1);
    mdl.B(1,:) = 1; % Set all rows to ones as a bad initial condition
    mdl.B(:, end) = 0; % Set the last column (constant x_offset) to zero
    mdl.C = zeros(1, mdl.order);
    mdl.C(1) = 1;
    mdl.D = zeros(1, mdl.num_inputs+1);
    if mdl.delay_per_chan
        mdl.delay = zeros(mdl.num_inputs, 1);
    end
end

function sys = makesys(mdl)    
    mdl.A(:, 1) = -(mdl.poles); % Set first column to poles
    rhs = diag(ones(mdl.order - 1,1));
    rhs(end+1, :) = 0;
    mdl.A = [-abs(mdl.poles) rhs];
    %delayterms = struct('delay', abs(mdl.delay ./ 1000), ...
    %                           'a',[],'b', mdl.B, 'c', [],'d', []);    
    %sys = delayss(mdl.A, zeros(size(mdl.B)), mdl.C, mdl.D, delayterms);
    sys = ss(mdl.A, mdl.B, mdl.C, mdl.D);
    tmp = abs(mdl.delay) / 1000;
    tmp(end+1) = 0; % No delay for the constant offset
    if (mdl.delay_per_chan)
        sys.InputDelay = tmp;
    else
        sys.InputDelay(:) = abs(mdl.delay) / 1000; % (milliseconds)
    end    
end

function x = do_pole_zeros(mdl, x, stack, xxx)
    
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
             u(:, end+1) = 1; % Append on a constant 1 
             
             warning off Control:analysis:LsimStartTime;
             tmp(:,s) = lsim(sys, u, t, x0);
             warning on Control:analysis:LsimStartTime;
         end
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = tmp + mdl.y_offset;
    end
end

function do_plot_pz_impulse_response(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    impulse(sys);   
    do_xlabel('Time [s]');
    do_ylabel('Impulse Response');
end

function do_plot_pz_step_response(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};    
    sys = makesys(mdls{1});
    step(sys);
    do_xlabel('Time [s]');
    do_ylabel('Step Response');
    
end

function do_plot_pz_bodemag_plot(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    bodemag(sys);
    do_xlabel('Freq [Hz]');
    do_ylabel('Freq Response');
end

end