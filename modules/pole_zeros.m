function m = pole_zeros(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @pole_zeros;
m.name = 'pole_zeros';
m.fn = @do_pole_zeros;
m.pretty_name = 'Pole/Zeros';
m.editable_fields = {'poles', 'zeros', 'delays', 'gains', ...
                     'y_offset', 'input', 'time', 'output'};  
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.n_inputs    = 1;
m.n_poles     = 2; 
m.n_zeros     = 1;
m.delays    = [5]; % Input delays in ms
m.gains     = [1];
m.y_offset = 0;
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_zplane;
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
m.plot_fns{4}.fn = @do_plot_zplane;
m.plot_fns{4}.pretty_name = 'ZPlane Plot';

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
      
    mdl.n_inputs = C;
    mdl.gains  = ones(1, mdl.n_inputs);
    mdl.delays = 5 + zeros(mdl.n_inputs, 1);
    
    % This ad-hoc initialization works tolerably for n_zeros < 5
    mdl.poles = repmat(-50 + -10*[1:mdl.n_poles], mdl.n_inputs, 1); 
    mdl.zeros = repmat(0   + -10*[1:mdl.n_zeros], mdl.n_inputs, 1);
    
end

function sys = makesys(mdl)   
    p = {};
    z = {};
    for ii = 1:mdl.n_inputs
        p{ii} = mdl.poles(ii, :);
        z{ii} = mdl.zeros(ii, :);
    end    
    sys = zpk(z, p, mdl.gains);
    sys.InputDelay = abs(mdl.delays) / 1000; % (milliseconds)
end

function x = do_pole_zeros(mdl, x, stack, xxx)    
    sys = makesys(mdl);    
    for sf = fieldnames(x.dat)', sf=sf{1};        
         [T, S, C] = size(x.dat.(sf).(mdl.input));         
         tmp = zeros(T, S, 1);         
         for s = 1:S            
             t = x.dat.(sf).(mdl.time)(:,1);
             u = squeeze(x.dat.(sf).(mdl.input)(:, s, :));         
             u_nan = isnan(u);
             u(u_nan) = 0; % TODO: REPLACE THIS HACK?
             u(isinf(u)) = 10^6;
             warning off Control:analysis:LsimStartTime;
             tmp(:,s) = lsim(sys, u, t);
             warning on Control:analysis:LsimStartTime;
             tmp(:,u_nan) = nan;
         end
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

function do_plot_zplane(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    %z = [];
    %for ii = 1:length(sys.z)
    %    z(:,ii) = sys.z{ii};
    %    p(:,ii) = sys.p{ii};
    %end
    zplane(mdls{1}.zeros', mdls{1}.poles');
    do_xlabel('Real Axis');
    do_ylabel('Imaginary Axis');
end

end