function fit05cperfile()

global STACK XXX

% Remove an
if strcmp(STACK{end}{1}.name, 'correlation')
    STACK = STACK(1:end-1);
    XXX = XXX(1:end-1);
end

% Then boost on each file individually
disp('NOW FITTING EACH FILECODE SEPARATELY');
[~,mseidx]=find_modules(STACK,'mean_squared_error');
split_stack(2,mseidx{1}-1);
fit05c;

%fit_split_simply(@fit05c, @split_by_filecode, @unify_respfiles); 

