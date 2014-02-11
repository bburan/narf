function sparsecheat()
% This method "cheats" by peeking at the VALIDATION set performance. While
% this is cheating since it introduces the possibility of undetectable
% overfitting, it may provide useful information about what level of
% sparsity should be applied to a particular batch of neurons. 

global META;

fitter = @fit10;

fitter();

sparsity_penalties = logspace(-7, -1, 20);
META.sparsity_penalty_results = [];

for ii = 1:length(sparsity_penalties) 
    fprintf('===> SPARSECHEAT %d/20\n', ii);
    META.sparsity_weight = sparsity_penalties(ii);
    fitter();
    [~,~, val_nmse] = META.perf_metric();
    META.sparsity_penalty_results(ii) = val_nmse;
end
