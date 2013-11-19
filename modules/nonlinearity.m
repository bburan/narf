function m = nonlinearity(args)
% Applies a general nonlinear function to the input. 
%
% The nonlinear function (NLFN) must be a function of two arguments:
%
%    NLFN(PHI, X) -> Y
%
% where
%    PHI is NLFN's parameter vector, such as polynomial coefficients
%    X is the matrix that the function works on
%    Y is the output of NLFN
%    
% Good choices for NLFN are 
%    polyval         (See matlab documentation)
%    nl_log          (See utils/ directory)
%    nl_root         (See utils/ directory)
%    nl_sigmoid      (See utils/ directory)
%    nl_exponential  (See utils/ directory)
%    nl_zerothresh   (See utils/ directory)
%
% You may of course define your own functions. 

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonlinearity;
m.name = 'nonlinearity';
m.fn = @do_nonlinearity;
m.pretty_name = 'Generic Nonlinearity';
m.editable_fields = {'input_stim', 'input_resp', 'time', 'output', 'nlfn', 'phi'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.nlfn = @polyval;
m.phi = [1 0];   % Default is to pass the signal through untouched

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_smooth_scatter_and_nonlinearity;
m.plot_fns = {};

m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Output Channels (All)';

m.plot_fns{2}.fn = @do_plot_single_default_output;
m.plot_fns{2}.pretty_name = 'Output Channel (Single)';

m.plot_fns{3}.fn = @do_plot_scatter_and_nonlinearity; 
m.plot_fns{3}.pretty_name = 'Stim/Resp Scatter';

m.plot_fns{4}.fn = @do_plot_smooth_scatter_and_nonlinearity; 
m.plot_fns{4}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{5}.fn = @do_plot_channels_as_heatmap; 
m.plot_fns{5}.pretty_name = 'Output Channels (Heatmap)';

%m.plot_fns{1}.fn = @(stack, xxx) do_plot_nonlinearity(stack, xxx, stack{end}.input_stim, @(x) stack{end}.nlfn(stack{end}.phi, x), false);
%m.plot_fns{1}.pretty_name = 'Nonlinearity';

%m.plot_fns{4}.fn = @(stack, xxx) do_plot_nonlinearity(stack, xxx, stack{end}.input_stim, @(x) stack{end}.nlfn(stack{end}.phi, x), true);
%m.plot_fns{4}.pretty_name = 'Nonlinearity + Histogram';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_nonlinearity(mdl, x, stack, xxx)

    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        y = zeros(T, S, C);
        
        % Note: we assume nlfn is vector-valued
        y = mdl.nlfn(mdl.phi, x.dat.(sf).(mdl.input_stim));      
        
        % TODO: Find a better solution than this hacky way of zeroing nans
        % so that optimization continue in the presence of singularities
        y(isnan(y)) = 0;
        y(isinf(y)) = 10^6;

        x.dat.(sf).(mdl.output) = y;
    end
end



function help_plot_nonlinearity(sel, mdls, xins, xouts)     
    xlims = xlim();  
    xs = linspace(xlims(1), xlims(2), 100);   
       
    for ii = 1:length(mdls)    
        ys = mdls{ii}.nlfn(mdls{ii}.phi, xs);         
        xouts{ii}.dat.(sel.stimfile).input = xs';
        xouts{ii}.dat.(sel.stimfile).output = ys';
    end
    
    hold on;
    do_plot(xouts, 'input', 'output', ...
            sel, 'Nonlinearity Input [-]', 'Prediction [Hz]');   
	hold off;
    
end

function do_plot_scatter_and_nonlinearity(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));    
    do_plot_scatter(sel, xins, mdls{1}.input_stim, mdls{1}.input_resp, 100);  
    help_plot_nonlinearity(sel, mdls, xins, xouts);     
end

function do_plot_smooth_scatter_and_nonlinearity(sel, stack, xxx)            
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));    
    do_plot_scatter(sel, xins, mdls{1}.input_stim, mdls{1}.input_resp, 100);  
    help_plot_nonlinearity(sel, mdls, xins, xouts); 
end

function do_plot_channels_as_heatmap(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));  
    ii=1;  % assuming just a single data set for now...
    h = imagesc(xouts{ii}.dat.(sel.stimfile).(mdls{1}.time)(:),...
                1:size(xouts{ii}.dat.(sel.stimfile).(mdls{1}.output),3),...
                squeeze(xouts{ii}.dat.(sel.stimfile).(mdls{1}.output)(:, sel.stim_idx, :))');
    do_xlabel('Time [s]');
    do_ylabel('Output [-]');
    
    set(gca,'YDir','normal');
    axis xy tight;
end

end