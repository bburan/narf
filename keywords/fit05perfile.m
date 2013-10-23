function fit05perfile()

global STACK XXX

disp('FIRST SINGLE FIT FOR ALL DATA (NO SPLITTING)');
fit05;

% Then boost on each file individually
disp('THEN FIT EACH FILE SEPARATELY');
fit_split_simply(@fit05, @split_by_filecode, @unify_respfiles); 

