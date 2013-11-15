function m = pole_zeros(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @pole_zeros;
m.name = 'pole_zeros';
m.fn = @do_pole_zeros;
m.pretty_name = 'Pole/Zeros';
m.editable_fields = {'order', 'poles_real', 'poles_img', 'zeros_real', 'zeros_img', ...
                     'delays', 'gains', 'y_offset', 'input', 'time', 'output'};  
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.n_inputs    = 1;
m.order       = 2;
m.poles_real  = [1]; 
m.poles_img   = [0];
m.zeros_real  = [0];
m.zeros_img   = [0];
m.delays      = [5]; % Input delays in ms
m.gains       = [1];
m.y_offset    = 0;
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
      
    mdl.n_inputs = C;
    mdl.gains  = 10 * ones(mdl.n_inputs, 1);
    mdl.delays = 5 + zeros(mdl.n_inputs, 1);
    
    if ~(mod(mdl.order,2)  == 0)
        error('Order must be an even number because I hacked this together.');
    end
    
    mdl.poles_real  = randn(mdl.n_inputs, mdl.order/2); 
    mdl.poles_img   = randn(mdl.n_inputs, mdl.order/2);
    mdl.zeros_real  = randn(mdl.n_inputs, mdl.order);
    mdl.zeros_img   = []; % randn(mdl.n_inputs, mdl.order/2);            
    
end

function sys = makesys(mdl)   
    p = {};
    z = {};
    for ii = 1:mdl.n_inputs
        for jj = 1:mdl.order/2
            % Create poles in conjugate pairs            
            p{ii}(2*jj-1) = -abs(mdl.poles_real(ii,jj)) + 1i * abs(mdl.poles_img(ii,jj));
            p{ii}(2*jj) = -abs(mdl.poles_real(ii,jj)) - 1i * abs(mdl.poles_img(ii,jj));
            

            % TODO: Create zeros in conjugate pairs too? Is that physical?
            % z{ii} = [mdl.zeros_real + i*abs(mdl.zeros_img)...
            %          mdl.zeros_real - i*abs(mdl.zeros_img)];                        
            % z{ii}(jj+1) = mdl.zeros_real(ii);             
            % TODO: Try this?
            %  z{ii}(1) = [0];   % ALWAYS put a zero at DC frequency
            z{ii}(2*jj-1) = mdl.zeros_real(ii, 2*jj-1);  
            z{ii}(2*jj) = mdl.zeros_real(ii, 2*jj);  
        end
    end    
    sys = zpk(z, p, mdl.gains');
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
             % u(:, end+1) = 1; % Append on a constant 1 
             
             warning off Control:analysis:LsimStartTime;
             tmp(:,s) = lsim(sys, u, t);
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