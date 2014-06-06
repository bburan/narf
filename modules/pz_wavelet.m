function m = pz_wavelet(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @pz_wavelet;
m.name = 'pz_wavelet';
m.fn = @do_pz_wavelet;
m.pretty_name = 'PZ Wavelet';
m.editable_fields = {'center_freq_khz', 'Q_factor', 'N_order', ...
                     'M_order', 'zeros', 'poles', ...
                     'delayms', 'gain', 'y_offset', 'align_peak', ...
                     'input', 'time', 'output'};  
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.poles       = [];
m.zeros       = []; % In kHz
m.N_order     = 4;  % aka, number of pairs of poles
m.M_order     = 1;  % aka, number of stacked zeros
m.center_freq_khz = 2;
m.Q_factor    = 3; % Actually the offset from 1/sqrt(2)
m.delayms     = 5; % Input delays in ms
m.gain        = 1;
m.y_offset    = 0;
m.align_peak  = false; % Align the peak of the impulse response to t=0?

m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.auto_plot = @do_plot_pz_bodemag_plot;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_single_default_output;
m.plot_fns{1}.pretty_name = 'Output vs Time';
m.plot_fns{2}.fn = @do_plot_pz_impulse_response;
m.plot_fns{2}.pretty_name = 'Impulse Response';
m.plot_fns{3}.fn = @do_plot_pz_step_response;
m.plot_fns{3}.pretty_name = 'Step Response';
m.plot_fns{4}.fn = @do_plot_pz_bodemag_plot;
m.plot_fns{4}.pretty_name = 'Bode Mag. Plot';
m.plot_fns{5}.fn = @do_plot_zplane;
m.plot_fns{5}.pretty_name = 'ZPlane Plot';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

% ------------------------------------------------------------------------
% INSTANCE METHODS

function [sys, t_align] = makesys(mdl)
    CF = abs(mdl.center_freq_khz) * 1000;
    Q  = 1/sqrt(2) + 0.00001 + abs(mdl.Q_factor); %.00001 is for numerical stability
    N  = ceil(abs(mdl.N_order));
    
    % Calculate useful parameterizations
    w_c = 2*pi*CF;    
    theta = asin(1-(1/(2*Q^2)));
    w_n = w_c / sqrt(1-(1/(2*Q^2))); % Natural frequency in rad/sec

    % Create a set of 2nd order poles that will be repeated N times
    b   = w_c * cos(theta/2 + pi/4);
    w_r = w_c * sin(theta/2 + pi/4); % Using w_c instead of w_n...why?        
    poles = [-b + 1i*w_r; -b - 1i*w_r];
    
    % Replicate the number of poles depending on the order of the filter
    poles = repmat(poles, N, 1);
    
    if isempty(mdl.zeros)
        % Create the APGF to have a gain of 1 at DC
        sys = zpk([], poles, 1);
        sys = zpk([], poles, 1/dcgain(sys));
    else
        z = repmat(mdl.zeros * 1000, mdl.M_order, 1);
        sys = zpk(z, poles, 1);
        % Gain is 0 at DC for OZGF, so we just try to keep it "near" the
        % same levels by normalizing by the ringing frequency w_r. Why does
        % this get slightly erronious for very low Q factors?
        sys = zpk(z, poles, (1/cos(theta)^N) / abs(evalfr(sys, w_r*1j)));         
    end

    % Correct for the time delay introduced by the wavelet
    if isfield(mdl, 'time_align') && mdl.time_align
        [imp, time] = impulse(sys);
        [~, idx] = max(abs(hilbert(imp)));
        t_align = time(idx);
    else
        t_align = 0;
    end
    
    sys.InputDelay = (abs(mdl.delayms) / 1000); % (milliseconds)
end

function x = do_pz_wavelet(mdl, x)    
    [sys, t_align] = makesys(mdl);    
    for sf = fieldnames(x.dat)', sf=sf{1};        
         [T, S, C] = size(x.dat.(sf).(mdl.input));         
         if C ~= 1
             error('There must be only one input channel');
         end
          
         tmp = zeros(T,S,C);
        
         for s = 1:S
            t = x.dat.(sf).(mdl.time)(:,1);
            u = squeeze(x.dat.(sf).(mdl.input)(:, s, :));
             
            % If there are NANs in the input, treat them like 0's when
            % simulating with lsim, then NAN out the output later.
            % This is kind of an ugly hack. 
            u_nan = isnan(u);
            u(u_nan) = 0;
            u(isinf(u)) = 10^6;
            warning off Control:analysis:LsimStartTime;
            warning off Control:analysis:LsimUndersampled;
             
            blah = lsim(sys, u, t);
            blah2 = hilbert(blah);
            tmp(:,s) = abs(blah2);
 
            warning on Control:analysis:LsimUndersampled;
            warning on Control:analysis:LsimStartTime;
            nanidxs = any(u_nan,2);
            tmp(nanidxs,s) = nan;                           
                        
            % Shift everything by t_align bins
            if isfield(mdl, 'time_align') && mdl.time_align
                idx = find(t < t_align, 1, 'last');
                tmp(:,s) = [tmp(idx:end,s)' tmp(end)*ones(1, idx-1)];
            end
         end
          
         % Apply the post-filter function
         x.dat.(sf).(mdl.output) = mdl.gain * abs(tmp) + mdl.y_offset;         
    end
end

function do_plot_pz_impulse_response(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    h = impulseplot(sys);   
    setoptions(h, 'Grid', 'on');      
    do_xlabel('Time [s]');
    do_ylabel('Impulse Response (Before Time Alignment)');
end

function do_plot_pz_step_response(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};    
    sys = makesys(mdls{1});
    stepplot(sys);
    setoptions(h, 'Grid', 'on');
    do_xlabel('Time [s]');
    do_ylabel('Step Response');    
end

function do_plot_pz_bodemag_plot(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    h = bodeplot(sys);
    setoptions(h,'FreqUnits','Hz','PhaseVisible','off', 'Grid', 'on');
    do_xlabel('Freq [Hz]');
    do_ylabel('Freq Response');
    axis tight;
end

function do_plot_zplane(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    zplane(sys.z{1}, sys.p{1});
    do_xlabel('Real Axis');
    do_ylabel('Imaginary Axis');
end

end
