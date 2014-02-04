function fit05perfile()

global STACK XXX

% Then boost on each file individually
disp('THEN FIT EACH FILE SEPARATELY');
fit_split_simply(@fit05, @split_by_filecode, @unify_respfiles); 
%fit_split_simply(@fit05s07, @split_by_filecode, @unify_respfiles); 

