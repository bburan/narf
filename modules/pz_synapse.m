function m = pz_synapse(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @pz_synapse;
m.name = 'pz_synapse';
m.fn = @do_pz_synapse;
m.pretty_name = 'PZ Synapse';
m.editable_fields = {'prefn', 'prephi', 'postfn', 'postphi', ...
                     'poles', 'zeros', 'delayms', 'gain', 'y_offset', ...
                     'input', 'time', 'output'};  
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.prefn  =  @polyval;
m.prephi = [1 1];
m.postfn  =  @nl_softzero;
m.postphi = [1 1 0 0];
m.poles       = [-50];
m.zeros       = [];
m.delayms     = 5; % Input delays in ms
m.gain        = 10;
m.y_offset    = 0;

m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_pz_impulse_response;
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

function sys = makesys(mdl)
    sys = zpk(mdl.zeros, mdl.poles, mdl.gain);
    sys.InputDelay = abs(mdl.delayms) / 1000; % (milliseconds)
end

function x = do_pz_synapse(mdl, x, stack, xxx)    
    sys = makesys(mdl);    
    for sf = fieldnames(x.dat)', sf=sf{1};        
         [T, S, C] = size(x.dat.(sf).(mdl.input));         
         tmp = zeros(T, S, 1);                                  
                  
         for s = 1:S            
             t = x.dat.(sf).(mdl.time)(:,1);
             u = squeeze(x.dat.(sf).(mdl.input)(:, s, :));
             
             % Apply a pre-filter function
             u = mdl.prefn(mdl.prephi, u);
             
             % If there are NANs in the input, treat them like 0's when
             % simulating with lsim, then NAN out the output later.
             % This is kind of an ugly hack. 
             u_nan = isnan(u);
             u(u_nan) = 0;
             u(isinf(u)) = 10^6;
             warning off Control:analysis:LsimStartTime;
             tmp(:,s) = lsim(sys, u, t);
             
             warning on Control:analysis:LsimStartTime;
             nanidxs = any(u_nan,2);
             tmp(nanidxs,s) = nan;
                          
         end
         % Apply the post-filter function
         x.dat.(sf).(mdl.output) = mdl.postfn(mdl.postphi, tmp) + mdl.y_offset;
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
    zplane(mdls{1}.zeros', mdls{1}.poles');
    do_xlabel('Real Axis');
    do_ylabel('Imaginary Axis');
end

end