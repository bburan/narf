function mixfit1()

mse();       % 
qboost();    % We think this 'sets up' initial conditions
boostirel(); % Boost on each module for efficiency and avoiding minima
qfmin();     % Another minima-resistant method
qlsqi();     % Zoom in on minima
qboosti();   % Slightly wider search
qboost();    % Final search
