function m = sparse_empirical_nonlinearity(args)
% Generates a sparse, smooth representation of a nonlinearity using a small
% number of gaussian mixtures.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @sparse_empirical_nonlinearity;
m.name = 'sparse_empirical_nonlinearity';
m.fn = @do_se_nonlinearity;
m.pretty_name = 'Sparse Empirical Nonlinearity';
m.editable_fields = {'relvar', 'input_stim', 'input_resp', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.get_parms = @init_nonparm_nonlinearity;
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.relvar = 0.25;  % Magnitude of the gaussian variance, relative to the 
                  % total range (max - min) of the input space.

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_smooth_scatter_and_nonlinearity; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{2}.fn = @(stack, xxx) do_plot_signal(stack, xxx, stack{end}.time, stack{end}.output);
m.plot_fns{2}.pretty_name = 'Output vs Time';

m.plot_fns{3}.fn = @do_plot_scatter_and_nonlinearity; 
m.plot_fns{3}.pretty_name = 'Stim/Resp Scatter';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function nl = calc_sparse_nonlinearity(stack, xxx)
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
        
    % Sort and average the entries so there are ~1000 data points
    D = [pred resp]; 
    D = sortrows(D);
    N = size(D, 1);
    D = conv_fn(D, 1, @nanmean, ceil(N/1000), 0);
    
    % TODO: This basis is incorrect because it assumes that the X are all
    % equally spaced apart, which is untrue
    function D2 = distSquared(X,Y)
        nx	= size(X,1);
        ny	= size(Y,1);
        D2 = (sum((X.^2), 2) * ones(1,ny)) + ...
             (ones(nx, 1) * sum((Y.^2),2)') - 2*X*Y';
    end
    
    X = D(:,1);
    Y = D(:,2);
    basisWidth	= (max(X) - min(X)) * mdl.relvar;
    basis = exp(-distSquared(X,X)/(basisWidth^2));
    [P, H, D] = SparseBayes('Gaussian', basis, Y);
    w = zeros(size(basis,2),1); % Basis weights
    w(P.Relevant) = P.Value;     
    Z = basis*w;
    Cx = X(P.Relevant);
    Cy = Y(P.Relevant);
    nl.Z = Z;
    nl.Cx = Cx;
    nl.Cy = Cy;
    nl.Q = H.Alpha;
    nl.X = X;
    nl.Y = Y;
    
    % We also return a function that can be used to interpolate arbitrary
    % values of y using the gaussian mixture
    function z = interpolomatic(xs)

        % TODO: Vectorize this function to do the gaussian mixture stuff
        
        %l = length(xs);
        %dists = abs(xs * ones(1, length(Cx)) - ones(l,1) * Cy');
        %z = sum(exp(-dists/(basisWidth^2)), 2);
        
        %dists = distSquared(xs, Cx);
        %z = sum(exp(-dists/(basisWidth^2)), 2);
                
        % VERY SLOW BUT SAFE WAY
        for ii = 1:length(xs)
            x = xs(ii);
            i1 = find(x > X, 1, 'last');   
            i2 = find(x <= X, 1, 'first');
            if isempty(i1)  % Below lower bound case
                z(ii) = Z(1);
            elseif isempty(i2) % Above upper bound case
                z(ii) = Z(end);
            else % Bounded on both sides case
                d = (x - X(i1)) / (X(i2) - X(i1));
                z(ii) = (1-d)*Z(i1) + d*Z(i2); 
            end
        end
    end
    nl.F = @interpolomatic;
 end

function x = do_se_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    nl = calc_sparse_nonlinearity(stack, xxx);
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        yout = zeros(T, S, C);
        yout = nl.F(x.dat.(sf).(mdl.input_stim)(:));
        x.dat.(sf).(mdl.output) = reshape(yout,[T,S,C]);
    end
end

function plotthecurve(nl)
    %hold on;
    %plot(nl.X, nl.Y, 'g.');
    plot(nl.X, nl.Z, 'b-');
    errorbar(nl.Cx, nl.Cy, nl.Q, 'r.');
    axis([min(nl.X), max(nl.X), min(nl.Y), max(nl.Y)]);
    %hold off;
end

function do_plot_scatter_and_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    nl = calc_sparse_nonlinearity(stack, xxx(1:end-1));
    hold on;
    do_plot_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
    plotthecurve(nl);
    hold off;
end

function do_plot_smooth_scatter_and_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    nl = calc_sparse_nonlinearity(stack, xxx(1:end-1));
    hold on;
    do_plot_avg_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);    
    plotthecurve(nl);
    hold off;
end

end