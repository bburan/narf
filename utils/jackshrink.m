function termcond = jackshrink(n_jacks, fitter, prior)
% Performs jackknifing by rewriting the XXX{2} data structure, which 
% is assumed to be produced by load_stim_from_baphy.
% ARGUMENTS:
%    N_JACKS    Number of jackknifes. Default: 10.
%    FITTER     Fitting routine to use. Default: fit_lsq
%    PRIOR      A prior belief to shrink towards. Default: all zeros.
% OPERATION:
%   Calls FITTER with XXX{2} set to each jackknife.
%   For the jackknife to work, FITTER must NOT alter XXX{2}. So any
%   recalc_xxx()'s that occur in FITTER must begin at 2 or higher.
%   Finally, it computes the appropriate shrinkage for the parameters using
%   the James Stein Estimator equation. 
%   http://upload.wikimedia.org/math/c/4/c/c4c103ca9fc386d90d9efeb64638f1f9.png
%   This allows it to work with ANY fitter and ANY performance metric.

global XXX STACK;

if nargin < 1
    n_jacks = 10;
end

if nargin < 2
    fitter = @fit_lsq;
end

recalc_xxx(1); 
cache = XXX{2};

phi_init = pack_fittables(STACK);


if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;
    return 
end

if nargin < 3
    prior = zeros(size(phi_init)); 
end

start_depth = find_fit_start_depth(STACK);

if start_depth < 2
    error('Regretfully, this implementation cannot jackknife parameters of STACK{1}.');
end

phi_jack = zeros(length(phi_init), n_jacks);

for jj = 1:n_jacks
    fprintf('\nJackknife [%d/%d]\n', jj, n_jacks);
    fns = fieldnames(XXX{2}.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        
        XXX{2}.dat.(sf) = cache.dat.(sf);
        
        % FIXME: Assumes only signals we need to jackknife stim, resp, and
        % respavg, which is not necessarily always true. 
        slen = size(cache.dat.(sf).stim, 1);
        rlen = size(cache.dat.(sf).resp, 1);
        ralen = size(cache.dat.(sf).respavg, 1);
        
        if ~isequal(slen, rlen) | ~isequal(rlen, ralen)
            error('Dim 1 lengths of stim, resp, respavg must be equal to jackknife.');
        end

        jackidx = floor(linspace(1, slen, n_jacks+1));
                
        XXX{2}.dat.(sf).stim(jackidx(jj):jackidx(jj+1),:,:) = NaN;
        XXX{2}.dat.(sf).resp(jackidx(jj):jackidx(jj+1),:,:) = NaN;
        XXX{2}.dat.(sf).respavg(jackidx(jj):jackidx(jj+1),:) = NaN;
    end

    unpack_fittables(phi_init);
    recalc_xxx(2); 
    fitter();
    phi_jack(:, jj) = pack_fittables(STACK);
end

% Compute how much to shrink each parameter 
m = length(phi_init);
n = n_jacks;
mu = mean(phi_jack, 2);
sigma_sq = var(phi_jack, [], 2);
% James-Stein Estimator that shrinks towards the prior.
if m < 3    
    fprintf('WARNING: A James Stein Estimator is only better if there are 3 or more params. Using simple mean.');
    phi_best = mu;
else
    phi_best = prior + (mu - prior) .* (  1 - (((m-2) * (sigma_sq / n)) / sum((mu - prior).^2))); 
end
% Stephen's shrinkage eqn is essentially the same, but without prior.
%phi_best = real(mu .* sqrt(1 - ((sqrt(sigma_sq) / sqrt(n)) ./ mu).^2));

% Recalc all the way through
unpack_fittables(phi_best);
XXX{2} = cache;
recalc_xxx(2);
termcond = NaN;

end
