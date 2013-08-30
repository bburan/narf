function mixfit1()

mse();

% alternative to qboost, run with quicker stop limit
global STACK;
phi_init = pack_fittables(STACK);
%fit_boost(n_steps, minstepsize, min_scoredelta, relative_delta)
fit_boost(length(phi_init),10^-9,10^-2,true);

boostirelsvd(); % Boost on each module for efficiency and avoiding minima
             %qfmin();     % Another minima-resistant method
             %qlsq();      % Zoom in on a minima
qboost();    % Final search of nearby minima
