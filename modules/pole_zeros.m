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
m.sum_channels = true; % When true, 
m.n_inputs    = 1;
m.n_poles     = 2; 
m.n_zeros     = 1;
m.delays      = [5]; % Input delays in ms
m.gains       = [1];
m.y_offset    = 0;
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.auto_plot = @do_plot_pz_heat_impulse_response;
m.auto_init = @auto_init_pz;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_pz_heat_impulse_response;
m.plot_fns{1}.pretty_name = 'Heat Impulse';
m.plot_fns{2}.fn = @do_plot_pz_impulse_response;
m.plot_fns{2}.pretty_name = 'Impulse Response';
m.plot_fns{3}.fn = @do_plot_pz_step_response;
m.plot_fns{3}.pretty_name = 'Step Response';
m.plot_fns{4}.fn = @do_plot_pz_bodemag_plot;
m.plot_fns{4}.pretty_name = 'Bode Mag. Plot';
m.plot_fns{5}.fn = @do_plot_zplane;
m.plot_fns{5}.pretty_name = 'ZPlane Plot';
m.plot_fns{6}.fn = @do_plot_single_default_output;
m.plot_fns{6}.pretty_name = 'Output vs Time';
% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

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
    mdl.gains  = 16 * ones(1, mdl.n_inputs);
    mdl.delays = 5 + zeros(mdl.n_inputs, 1);
    
    % This ad-hoc initialization works tolerably for n_zeros < 5
    mdl.poles = repmat(-30 + -20*[1:mdl.n_poles], mdl.n_inputs, 1); 
    mdl.zeros = repmat(-10 + -10*[1:mdl.n_zeros], mdl.n_inputs, 1);
    
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

function x = do_pole_zeros(mdl, x)    
    sys = makesys(mdl);    
    for sf = fieldnames(x.dat)', sf=sf{1};        
         [T, S, C] = size(x.dat.(sf).(mdl.input));     
         
         % Optimized to be 5x faster than previously.
         tmp = zeros(T*S, 1);         
         u = reshape(x.dat.(sf).(mdl.input), T*S, C);
         dt = x.dat.(sf).(mdl.time)(2) - x.dat.(sf).(mdl.time)(1);
         t = dt * [0:size(u,1)-1];
         u_nan = isnan(u);
         u(u_nan) = 0;
         u(isinf(u)) = 10^6;
         % warning off Control:analysis:LsimStartTime;
         tmp = lsim(sys, u, t);
         % warning on Control:analysis:LsimStartTime;
         tmp(u_nan) = nan;
         
         if prod(size(tmp))~=T*S,
             warning('pole_zeros: output size mismatch!?!?!?');
             tmp=tmp(1:size(u,1),:);
         end
         tmp = reshape(tmp, T, S, 1);
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
    zplane(mdls{1}.zeros(:)', mdls{1}.poles(:)');
    do_xlabel('Real Axis');
    do_ylabel('Imaginary Axis');
end

function do_plot_pz_heat_impulse_response(sel, stack, xxx)
    mdls = stack{end};    
    xins = {xxx(1:end-1)};        
    sys = makesys(mdls{1});
    
    x = xins{1};
    x = x{end};
    fns = fieldnames(x.dat);
    time = x.dat.(fns{1}).(mdls{1}.time);
    dt = time(2) - time(1);
    Y = impulse(sys, linspace(0, 0.15, 500));  
    Y = squeeze(Y);
    
    [mod, mod_idx] = find_modules(stack, 'weight_channels', true);
    
    if ~isempty(mod)
        img = Y' ./ mean(abs(Y(:)));
        img = cat(1, nan(1, size(img,2)), img);      
        [~, weights] = mod{1}.fn(mod{1}, xins{1}{mod_idx});
        if size(weights, 2) == size(Y,2)
            bot = weights * Y';
            img = cat(1, bot ./ mean(abs(bot(:))), img);        
            h = imagesc(img);
            do_ylabel('STRF | Impulse');    
        else
            h = imagesc(Y);    
            do_ylabel('Impulse Response');
        end
    else
        h = imagesc(Y);    
        do_ylabel('Impulse Response');
    end
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    set(gca,'YDir','normal');
    axis xy tight;   
    do_xlabel('Time [ms]');
   
    setAxisLabelCallback('X', @(z) z*150/500);
        
    
    
    
end


end
