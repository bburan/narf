function fit05cperfile()

global STACK XXX

% Then boost on each file individually
disp('NOW FITTING EACH FILECODE SEPARATELY');
fit_split_simply(@fit05c, @split_by_filecode, @unify_respfiles); 

