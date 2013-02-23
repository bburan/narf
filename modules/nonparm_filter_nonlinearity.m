function m = nonparm_filter_nonlinearity(args)
% Generates a quick, nonparametric nonlinearity by filtering the output
% with a wide gaussian and extrapolating beyond those points. 

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonparm_filter_nonlinearity;
m.name = 'nonparm_filter_nonlinearity';
m.fn = @do_npfnl;
m.pretty_name = 'Nonparametric Filter Nonlinearity';
m.editable_fields = {'gwinval', 'input_stim', 'input_resp', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.gwinval = 5;  % Good values probably between 3 and 6. Determines gaussian filter window size relative to total number of data point samples.
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_smooth_scatter_npfnl; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{2}.fn = @(stack, xxx) do_plot_signal(stack, xxx, stack{end}.time, stack{end}.output);
m.plot_fns{2}.pretty_name = 'Output vs Time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function npfnl = calc_npfnl(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Create the pred and resp vectors
    pred = []; 
    resp = [];
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};
        pred = cat(1,pred,x.dat.(sf).(mdl.input_stim)(:));
        resp = cat(1,resp,x.dat.(sf).(mdl.input_resp)(:));
    end
    
    % Remove entries with NaNs in resp.
    keepidx=find(~isnan(resp));
    pred=pred(keepidx);
    resp=resp(keepidx);
    
    % Build the gaussian smoothing filter   
    winsize = ceil(length(pred)/mdl.gwinval);
    if mod(winsize, 2) == 1
        winsize = winsize + 1;
    end
    gf = gausswin(winsize);
    gf = gf / sum(gf);
    
    % Throw out duplicates, sort, and smooth
    [~, idxs, ~] = unique(pred, 'first'); % TODO: Take mean instead of just throwing duplicates away
    D = sortrows([pred(idxs) resp(idxs)]); 
    S = iconv(D(:, 1), gf);
    R = iconv(D(:, 2), gf);
    
    right_slope = (R(end) - R(end-1)) / (S(end) - S(end-1));
    right_interp = @(x) R(end) + right_slope * (x-S(end));
    left_slope = (R(2) - R(1)) / (S(2) - S(1));
    left_interp = @(x) R(1) - left_slope * (S(1)-x);
    
    % Return a function that can interpolate arbitrary values
    function z = hacky_interpolator(xs)        
        
        z = interp1(S, R, xs, 'linear');
        
        % Correct the fact that values outside the linear range are NaN
        leftz  = xs < S(1); 
        rightz = xs > S(end);

        z(leftz) = left_interp(xs(leftz));
        z(rightz) = right_interp(xs(rightz));
         
    end
    npfnl = @hacky_interpolator;
 end

function x = do_npfnl(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    npfnl = calc_npfnl(stack, xxx);
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        yout = npfnl(x.dat.(sf).(mdl.input_stim)(:));
        x.dat.(sf).(mdl.output) = reshape(yout,[T,S,C]);
    end
end

function plotthecurve(nlfn, xs)
    xmin = min(xs);
    xmax = max(xs);
    aa = linspace(xmin, xmax, 100);
    z = nlfn(aa);
    ymin = min(z);
    ymax = max(z);
    plot(aa, z, 'k-');  
    axis([xmin, xmax, ymin, ymax]);
end

function do_plot_smooth_scatter_npfnl(stack, xxx)
    mdl = stack{end};
    
    npfnl = calc_npfnl(stack, xxx(1:end-1));
    
    hold on;
    do_plot_avg_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);    
    
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    plotthecurve(npfnl, xxx{end-1}.dat.(sf).(mdl.input_stim)(:));
    do_plot_avg_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp); 
    
    hold off;
end

end