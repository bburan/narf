function mixfit3perfile()

mse(); %% We need a performance metric, so add one
mixfit3nomse(); %% Relative boost algorithm across all files to initialize

% Then boost on each file individually
fit_split_simply(@mixfit3nomse, @split_by_respfile, @unify_respfiles); 
