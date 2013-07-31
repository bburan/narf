function mixfit0()

mse();       % 
qboost();    % We think this 'sets up' initial conditions
boostirel(); % Boost on each module for efficiency and avoiding minima
qfmin();     % Another minima-resistant method
qlsq();      % Zoom in on a minima
qboost();    % Final search of nearby minima