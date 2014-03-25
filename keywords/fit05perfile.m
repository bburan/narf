function fit05perfile()

global STACK XXX

% Then boost on each file individually
disp('NOW FITTING EACH FILECODE SEPARATELY');
fit_split_simply(@fit05, @split_by_filecode, @unify_respfiles); 

