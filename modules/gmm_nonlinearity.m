function m = gmm_nonlinearity(args)
% Makes a smooth nonlinearity using a gaussian mixture model.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @gmm_nonlinearity;
m.name = 'gmm_nonlinearity';
m.fn = @do_gmm_nonlinearity;
m.pretty_name = 'Gaussian Mixture Model';
m.editable_fields = {'num_pts', 'num_gaussians', 'input_stim', 'input_resp', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time       = 'stim_time';
m.output     = 'stim';
m.num_gaussians = 4; 
m.num_pts    = 100;
m.thresh     = 1e-10; 

% Optional fields
m.auto_init = @auto_init_gmm;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_smooth_scatter_and_gmm; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{2}.fn = @do_plot_scatter_and_gmm; 
m.plot_fns{2}.pretty_name = 'Stim/Resp Scatter';

m.plot_fns{3}.fn = @(stack, xxx) do_plot_signal(stack, xxx, stack{end}.time, stack{end}.output);
m.plot_fns{3}.pretty_name = 'Output vs Time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input_stim, m.input_resp, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

function Data = calc_data(mdl, x)
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
    S=pred(keepidx);
    R=resp(keepidx);
    % Remove duplicate stimulus values, if any
    % TODO: I should probably average duplicate values, not discard them
    [~, idxs, ~] = unique(S, 'first');
    
    % Sort the rows
    Data = sortrows([S(idxs), R(idxs)]);
    N = size(Data, 1);
    Data = conv_fn(Data, 1, @nanmean, ceil(N/mdl.num_pts), 0);
    Data = Data';
end

function mm = auto_init_gmm(stack, xxx)
    mm = m;
    
    Data = calc_data(mm, xxx{end});
        
    % Initialize gaussian mixture model starting point
    [Priors, Mu, Sigma] = EM_init_kmeans(Data, mm.num_gaussians);    
    
	mm.gmm_priors = Priors;
    mm.gmm_mu = Mu; 
    mm.gmm_sigma = Sigma;
end

function interpfn = calc_gmm_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    Data = calc_data(mdl, x);
    
    %disp(mdl.gmm_mu);
    
    tic;   
    if isfield(mdl, 'gmm_priors') && isfield(mdl, 'gmm_mu') &&...
       isfield(mdl, 'gmm_sigma') && isfield(mdl, 'thresh')
        [Priors, Mu, Sigma] = EM(Data, mdl.gmm_priors, mdl.gmm_mu, ...
                                mdl.gmm_sigma, mdl.thresh);
    else 
        [Priors, Mu, Sigma] = EM_init_kmeans(Data, mdl.num_gaussians);   
        [Priors, Mu, Sigma] = EM(Data, Priors, Mu, Sigma);
    end
    toc
    
    % Return a function that can be used to interpolate arbitrary values of
    % x using a gaussian mixture model.
    function [z, sig] = interpolomatic(xs)
        [z, sig] = GMR(Priors, Mu, Sigma, xs', [1], [2]);
    end
    
    % TODO: REMOVE THIS VERY NAUGHTY CODE! Updating the STACK in this way
    % is really frowned upon. But I can't think of another way of doing
    % this type of functionality. 
    global STACK;
    idx = length(stack);
   
    STACK{idx}.gmm_priors = Priors;
    STACK{idx}.gmm_mu = Mu;
    STACK{idx}.gmm_sigma = Sigma;
    
    interpfn = @interpolomatic;
 end

function x = do_gmm_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    nlfn = calc_gmm_nonlinearity(stack, xxx);

    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        yout = zeros(T, S, C);
        yout = nlfn(x.dat.(sf).(mdl.input_stim)(:));
        x.dat.(sf).(mdl.output) = reshape(yout,[T,S,C]);
    end
end

function plotthecurve(nlfn, xs)
    xmin = min(xs);
    xmax = max(xs);
    aa = linspace(xmin, xmax, 100);
    [z, sig] = nlfn(aa');
    ymin = min(z);
    ymax = max(z);
    %D = [aa, z];
    %plotGMM(D', sig(1,1,:), [0 0 .8], 3);
    plot(aa, z, 'k-');  
    axis([xmin, xmax, ymin, ymax]);
end

function do_plot_scatter_and_gmm(stack, xxx)
    mdl = stack{end};
    
    nlfn = calc_gmm_nonlinearity(stack, xxx(1:end-1));
    hold on;
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    plotthecurve(nlfn, xxx{end-1}.dat.(sf).(mdl.input_stim)(:));    
    do_plot_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
    hold off;
end

function do_plot_smooth_scatter_and_gmm(stack, xxx)
    mdl = stack{end};
    
    nlfn = calc_gmm_nonlinearity(stack, xxx(1:end-1));
    hold on;
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    do_plot_avg_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
    plotthecurve(nlfn, xxx{end-1}.dat.(sf).(mdl.input_stim)(:));    
    hold off;
end

end