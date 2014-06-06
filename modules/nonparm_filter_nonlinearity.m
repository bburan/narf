function m = nonparm_filter_nonlinearity(args)
% Generates a quick, nonparametric nonlinearity by filtering the output
% with a wide gaussian and extrapolating beyond those points. 

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonparm_filter_nonlinearity;
m.name = 'nonparm_filter_nonlinearity';
m.fn = @do_npfnl;
m.pretty_name = 'Nonparametric Filter Nonlinearity';
m.editable_fields = {'gwinval', 'leftzero', 'edgesize', 'input_stim', 'input_resp', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.gwinval = 4;  % Good values probably between 3 and 6. 
                % Determines gaussian filter window size relative to total number of data point samples.
m.leftzero = false;
m.edgesize = 20; % Arbitrary number. Usually about 20 is good.
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
    
% TODO: These are magic numbers! 20 points and 80 points just seemed to
% work pretty well, but its' really ad-hoc!
    
% Optional fields
m.auto_plot = @do_plot_smooth_scatter_npfnl;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_smooth_scatter_npfnl; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{2}.fn = @do_plot_scatter_npfnl; 
m.plot_fns{2}.pretty_name = 'Stim/Resp Scatter';

m.plot_fns{3}.fn = @do_plot_single_default_output;
m.plot_fns{3}.pretty_name = 'Output vs Time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input_stim, m.input_resp, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

function npfnl = calc_npfnl(mdl, x)   
    % Create the pred and resp vectors
    pred = []; 
    resp = [];
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};
        pred = cat(1, pred, x.dat.(sf).(mdl.input_stim)(:));
        resp = cat(1, resp, x.dat.(sf).(mdl.input_resp)(:));
    end
    
    if (length(pred) ~= length(resp))
        error('Length of pred and resp must be equal!');
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
    D = excise(D);
    S = iconv(D(:, 1), gf);
    R = iconv(D(:, 2), gf);
    
    rs = mdl.edgesize; 
    ls = mdl.edgesize * 4; 
    
    if length(S) < mdl.edgesize * 4 || length(R) < mdl.edgesize * 4
        right_slope = 0;
        left_slope = 0;
    else
        right_slope = (mean(R(end-rs:end)) - mean(R(end-ls:end-rs))) / ...
                      (mean(S(end-rs:end)) - mean(S(end-ls:end-rs)));
        left_slope = (mean(R(1:rs)) - mean(R(rs+1:ls))) / ...
                     (mean(S(1:rs)) - mean(S(rs+1:ls)));
    end
    
	if mdl.leftzero
        left_slope = 0;
    end
    
    right_interp = @(x) R(end) + right_slope * (x-S(end));
    left_interp = @(x) R(1) - left_slope * (S(1)-x);
    
    lengthsok = (length(S) > 1) && (length(R) > 1);
    
    % Return a function that can interpolate arbitrary values
    function z = hacky_interpolator(xs)        
        
        if lengthsok
            % Quick fix is just to catch the error
            try
                z = interp1(S, R, xs, 'linear');
            catch err
                z = xs;
            end            
        else
            z = xs;
        end
        
        % Correct the fact that values outside the linear range are NaN
        leftz  = xs < S(1); 
        rightz = xs > S(end);

        z(leftz) = left_interp(xs(leftz));
        z(rightz) = right_interp(xs(rightz));
         
    end
    npfnl = @hacky_interpolator;
 end

function x = do_npfnl(mdl, x)
    
    npfnl = calc_npfnl(mdl, x);
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        in = x.dat.(sf).(mdl.input_stim);
        [T, S, C] = size(in);
        if all(isnan(in(:)))
            yout = in(:);
        else
            yout = npfnl(in(:));
        end
        x.dat.(sf).(mdl.output) = reshape(yout,[T,S,C]);
    end
end

function help_plot_npfnl(sel, mdls, xins, xouts)     
    
    v = axis;
    xmin = v(1);
    xmax = v(2);
    
    for ii = 1:length(mdls)    
        npfnl = calc_npfnl(mdls{ii}, xins{ii}{end});
        
        xs = xins{ii}{end}.dat.(sel.stimfile).(mdls{ii}.input_stim)(:);  
        xmin = min(xs);
        xmax = max(xs);
        aa = linspace(xmin, xmax, 100);        
        z = npfnl(aa);
        
        xouts{ii}.dat.(sel.stimfile).forprinting_in = aa;
        xouts{ii}.dat.(sel.stimfile).forprinting_out = z';
    end
    
    hold on;
    do_plot(xouts, 'forprinting_in', 'forprinting_out', ...
            sel, 'NPFNL Input [-]', 'RespAvg Prediction [Hz]');   
	%hold off;
end

function do_plot_smooth_scatter_npfnl(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end));    
    %do_plot_scatter(sel, xins, mdls{1}.input_stim, mdls{1}.input_resp, 100);  
    help_plot_npfnl(sel, mdls, xins, xouts); 
end

function do_plot_scatter_npfnl(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end));    
    do_plot_scatter(sel, xins, mdls{1}.input_stim, mdls{1}.input_resp);  
    help_plot_npfnl(sel, mdls, xins, xouts);
end

end
