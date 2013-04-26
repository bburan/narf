function score = pm_mse()
% pm_mse()
%
% Performance metric: MSE with no sparseness penalties

global XXX;

score = XXX{end}.score_train_mse; % TODO: Assumes existence of MSE module