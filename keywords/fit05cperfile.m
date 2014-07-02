function fit05cperfile()

global STACK XXX

disp('NOW FITTING EACH TRIAL_CODE SEPARATELY');

% Remove any correlation at end of stack
if strcmp(STACK{end}{1}.name, 'correlation')
    STACK = STACK(1:end-1);
    XXX = XXX(1:end-1);
end

% find first STACK entry with fit_fields
ii=1;
while ii<length(STACK) && ...
        (~isfield(STACK{ii}{1},'fit_fields') || ...
         isempty(STACK{ii}{1}.fit_fields)),
    ii=ii+1;
end

[~,mseidx]=find_modules(STACK,'mean_squared_error');

% Then boost on each file individually
split_stack(ii,mseidx{1}-1);
fit05c;


