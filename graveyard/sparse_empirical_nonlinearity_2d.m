function m = sparse_empirical_nonlinearity_2d(args)
% Generates a sparse, smooth representation of a nonlinearity using a small
% number of gaussian mixtures.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @sparse_empirical_nonlinearity_2d;
m.name = 'sparse_empirical_nonlinearity_2d';
m.fn = @do_se_nonlinearity_2d;
m.pretty_name = 'Sparse Empirical Nonlinearity 2D';
m.editable_fields = {'numpts', 'relvar', 'input_stim', 'input_resp', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.numpts = 100; % Doubling number of points takes 4x the calc time
                % However, more points is usually better
                % Good range seems to be 50 to 1000 points.
                % 100 is pretty fast and 'good enough'
m.relvar = 0.2;  % Magnitude of the gaussian variance, relative to the 
                 % total range (max - min) of the input space.

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_smooth_scatter_and_nonlinearity_2d; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{2}.fn = @(stack, xxx) do_plot_signal(stack, xxx, stack{end}.time, stack{end}.output);
m.plot_fns{2}.pretty_name = 'Output vs Time';

m.plot_fns{3}.fn = @do_plot_scatter_and_nonlinearity_2d; 
m.plot_fns{3}.pretty_name = 'Stim/Resp Scatter';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function nl = calc_sparse_nonlinearity_2d(stack, xxx)
    % Returns:
    %   X : X-axis points, taken emperically. 
    %   Z : Predicted curve at each value in X
    %   C : Gaussian control points (centers)
    %   Q : Precision of the control points (inv of variance)
    mdl = stack{end};
    x = xxx{end};
    nl = [];
    % Create the pred and resp vectors
    pred = []; 
    resp = [];
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};
        pred = cat(1,pred,x.dat.(sf).(mdl.input_stim)(:));
        resp = cat(1,resp,x.dat.(sf).(mdl.input_resp)(:));
    end
    
    % Remove entries with NaNs in resp.
    % TODO: If we want to do shrinkage with this fitter, we would also need
    % to remove NaNs from the stimulus. 
    keepidx=find(~isnan(resp));
    pred=pred(keepidx);
    resp=resp(keepidx);
        
    % Sort and average the entries so there are data points
    D = [pred resp]; 
    D = sortrows(D);
    N = size(D, 1);
    D = conv_fn(D, 1, @nanmean, ceil(N/mdl.numpts), 0);
    
    X = D(:,1);
    Y = D(:,2);
    XbasisWidth	= (max(X) - min(X)) * mdl.relvar;
    YbasisWidth	= (max(Y) - min(Y)) * mdl.relvar;
    
    function D2 = distSquared(Xi,Yi)
        nx	= size(Xi,1);
        ny	= size(Yi,1);
        D2 = (sum((Xi.^2), 2) * ones(1,ny)) + ...
             (ones(nx, 1) * sum((Yi.^2),2)') - 2*Xi*Yi';
    end    

    Xbasis = exp(-distSquared(X,X)/(XbasisWidth^2));
    Ybasis = exp(-distSquared(Y,Y)/(YbasisWidth^2));
    [PX, HX, DX] = SparseBayes('Gaussian', basis, X);
    w = zeros(size(Xbasis,2),1); % Basis weights
    w(PX.Relevant) = PX.Value;     
    ZX = Xbasis*w;
    CxX = X(PX.Relevant);
    CyX = Y(PX.Relevant);
    nl.Z = ZX;
    nl.CxX = CxX;
    nl.CyX = CyX;
    nl.QX = HX.Alpha;
    nl.X = X;
    nl.Y = Y;
    
%     % Compute an extra control point to the right
%     slope    = (P.Value(end) - P.Value(end-1)) / (Cx(end) - Cx(end-1));
%     Cx_extra = Cx(end) + (Cx(end) - Cx(end-1));
%     Cw_extra = P.Value(end) + (slope * (Cx_extra - Cx(end)));
    
    % We also return a function that can be used to interpolate arbitrary
    % values of y using the gaussian mixture
    function z = interpolomatic_2d(xs)
        % TODO: Vectorize this function to do the gaussian mixture stuff
        newbasis = exp(-distSquared(xs, CxX) / (XbasisWidth^2));
        z = newbasis * PX.Value;
        
        % If we don't interpolate beyond the bounds of the training set,
        % they quickly fall off to zero. These 'stretch' the starting and
        % ending values, which is a very stupid approximation but also is
        % pretty safe. 
        % tops = xs > max(X);
        % slope = (Cy(end) - Cy(end-1)) / (Cx(end) - Cx(end-1));
        % start = Cy(end);
        % z(tops) = xs(tops) + Cy(end) + ;
         
        % On the bottom side, it's simplest just to let the values fall of
        % to zero as they go sufficiently below any observed training set
        % values. We don't have to do anything for that to happen, gaussian
        % mixture models have this property already. 

    end
    nl.F = @interpolomatic_2d;
 end

function x = do_se_nonlinearity_2d(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    nl = calc_sparse_nonlinearity_2d(stack, xxx);
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        yout = zeros(T, S, C);
        yout = nl.F(x.dat.(sf).(mdl.input_stim)(:));
        x.dat.(sf).(mdl.output) = reshape(yout,[T,S,C]);
    end
end

function plotthecurve_2d(nl, xmin, xmax)
    
    %plot(nl.X, nl.Y, 'g.');
    %plot(nl.X, nl.Z, 'r-');
    %plot(nl.X, nl.F(nl.X), 'm-');
    aa = linspace(xmin, xmax, 100);
    plot(aa', nl.F(aa'), 'b-');
    plot(nl.Cx, nl.Cy, 'rx');
    % errorbar(nl.Cx, nl.Cy, nl.Q, 'r.');
    % axis([min(nl.X), max(nl.X), min(nl.Y), max(nl.Y)]);
    axis tight;

end

function do_plot_scatter_and_nonlinearity_2d(stack, xxx)
    mdl = stack{end};
    
    nl = calc_sparse_nonlinearity_2d(stack, xxx(1:end-1));
    hold on;
    do_plot_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
    
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    xmin = min(xxx{end-1}.dat.(sf).(mdl.input_stim)(:));
    xmax = max(xxx{end-1}.dat.(sf).(mdl.input_stim)(:));
    plotthecurve_2d(nl, xmin, xmax);
    
    hold off;
end

function do_plot_smooth_scatter_and_nonlinearity_2d(stack, xxx)
    mdl = stack{end};
    
    nl = calc_sparse_nonlinearity_2d(stack, xxx(1:end-1));
    hold on;
    do_plot_avg_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);    
    
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    xmin = min(xxx{end-1}.dat.(sf).(mdl.input_stim)(:));
    xmax = max(xxx{end-1}.dat.(sf).(mdl.input_stim)(:));    
    plotthecurve_2d(nl, xmin, xmax);
    hold off;
end

end