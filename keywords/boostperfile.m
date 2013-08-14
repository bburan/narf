function boostperfile()

mse(); %% We need a performance metric, so add one
booster(); %% Relative boost algorithm across all files to initialize

% Then boost on each file individually
fit_split_simply(@booster, @split_by_respfile, @unify_respfiles); 

