function fit05cperfile()

global STACK XXX

% Then boost on each file individually
disp('NOW FITTING EACH FILECODE SEPARATELY');
[~,mseidx]=find_modules(STACK,'mean_squared_error');
split_stack(2,mseidx{1}-1);
update_xxx(2);
fit05c;

%fit_split_simply(@fit05c, @split_by_filecode, @unify_respfiles); 

