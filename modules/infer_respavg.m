function m = infer_respavg(args)
% Infers the respavg

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @infer_respavg;
m.name = 'infer_respavg';
m.fn = @do_infer_respavg;
m.pretty_name = 'Infer Respavg';
m.editable_fields = {'input', 'time', 'output', 'strength'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'respavg'; 
m.time  = 'resp_time';
m.output = 'respavg';
m.distribution = 'poissexp'; % Options: normid, binomlogit, possexp
m.strength = 0.1;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_inferred_respavg;
m.plot_fns{1}.pretty_name = 'Training Set Respavg';
m.plot_gui_create_fn = @create_chan_selector_gui;

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

function x = do_infer_respavg(mdl, x)    
    glmopts.family = mdl.distribution;
    glmopts.algo = 'medium'; % 'medium', 'large', or 'newton'

    % Smooth ONLY the training set RESPAVG signals
    for ii = 1:length(x.training_set)
        sf = x.training_set{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        out = zeros([T, S, C]);
        
       %% TODO: Use RESP for the inference and then create a new RESPAVG 
       %% 
        
        for c = 1:C,
            for s = 1:S
                y = x.dat.(sf).(mdl.input)(:,s,c);
                y = ceil(y); % FIXME: This is just a hack for now        
                X = speye(length(y));
                P = spdiags(ones(size(y))*[-1,1], [0,1], length(y)-1, length(y));
                Q = P'*P;
                res = glmfitqp(y, X, Q*(mdl.strength), glmopts);
                out(:,s,c) = res.w;
            end
        end
        x.dat.(sf).(mdl.output) = out;
    end
end

function do_plot_inferred_respavg(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    
    hold on;
    do_plot(xouts, mdls{1}.time, mdls{1}.input, ...
            sel, 'Time [s]', 'Respavg [-]');
    hold off;     
end

end
