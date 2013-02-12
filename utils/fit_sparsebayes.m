function termcond = fit_sparsebayes(pred_field, target_field)

global XXX STACK;

if nargin < 2
    pred_field = 'stim';
end
if nargin < 2
    target_field = 'respavg';
end

start_depth = find_fit_start_depth(STACK);
phi_init = pack_fittables(STACK);
target = flatten_field(XXX{end}.dat, XXX{1}.training_set, target_field);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;
    return 
end

function basis = build_basis(phi)
    M = length(phi);
    N = length(target);    
       
    % Get the default prediction from where you are    
    recalc_xxx(start_depth);
    pred_default = flatten_field(XXX{end}.dat, XXX{1}.training_set, pred_field);
    
    % BASIS: NxM. Would be the prediction of each FIR coef (M coefs)?
    %        Each column is a basis vector
    % TARGETS: The output vector

    basis = zeros(N, M);
    
    % Perturb each fittable field by delta
    delta = 1;
    
    % Create basis vectors by subtracting default from the new prediction
    % The TARGET is respavg
    for jj = 1:length(phi)
        newphi = phi;
        newphi(jj) = phi(jj) + delta;
        unpack_fittables(newphi);
        recalc_xxx(start_depth);
        pred = flatten_field(XXX{end}.dat, XXX{1}.training_set, pred_field);
        pred(isnan(pred)) = 0;
        basis(:, jj) = pred - pred_default;
    end
    
end

fprintf('Fitting %d variables with SparseBayes()\n', length(phi_init));
N_iterations = 500;    

for ii = 1:N_iterations
    fprintf('SparseBayes Iteration %d/%d\n', ii, N_iterations);
    phi = pack_fittables(STACK);
    basis = build_basis(phi);
    [P, H, D] = SparseBayes('Gaussian', basis, target);
    phi_delta = zeros(size(phi));
    phi_delta(P.Relevant) = P.Value;    
    unpack_fittables(phi + phi_delta);
end

termcond = NaN;

end