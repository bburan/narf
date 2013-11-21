function m = nim_wrapper(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nim_wrapper;
m.name = 'nim_wrapper';
m.fn = @do_nim_wrapper;
m.pretty_name = 'NIM Wrapper';
m.editable_fields = {'input_stim', 'input_resp', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.nim_model = [];

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_nim_filters;
m.plot_fns{1}.pretty_name = 'NIM Filters';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_nim_wrapper(mdl, x, stack, xxx)   
    
    NMMPATH = '/auto/user/ivar/matlab/nmm/';
    addpath(genpath(NMMPATH));
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        % TODO: Presently, this only supports a single channel
        stim = x.dat.(sf).(mdl.input_stim)(:,:,1);
        stim = stim(:);        
        meanstd(1) = sqrt(mean(stim.^2));
        stim(:,1) = stim(:,1)/meanstd(1);        
        SR = 200; 
        tent_basis_spacing = 1;
        nLags = 30; 
        frac = 2;
        params_stim = NMMcreate_stim_params( nLags, 1/SR, frac, tent_basis_spacing );
        Xstim = create_time_embedding(stim, params_stim);

        [LL, penLL, pred_rate]= NMMCmodel_eval(mdl.nim_model, [], Xstim);
        
        x.dat.(sf).(mdl.output) = reshape(decimate(pred_rate,2) * SR, T,S);
       
    end
end

function do_plot_nim_filters(sel, stack, xxx)
    mdls = stack{end};
    xins = {xxx(1:end-1)};
    mdl = mdls{1};
    x = xins{1};
    
    hold on;
    h1= stem(mdl.nim_model.mods(1).filtK, 'g');
    h2= stem(mdl.nim_model.mods(2).filtK, 'b'); 
    hold off;        
    axis tight;
    do_xlabel('Coef Time Index');
    do_ylabel('Coef Magnitude');
    legend({'Excitatory', 'Inhibitory'}, 'Location', 'NorthWest');
end

end