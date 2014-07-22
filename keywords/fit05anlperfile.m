% function fit05anlperfile()
%
% use SEMSE rather than standard MSE
function fit05anlperfile()

global STACK XXX

disp('NOW FITTING POST-FIR STAGES PER FILE');

% Remove any correlation at end of stack
if strcmp(STACK{end}{1}.name, 'correlation')
    STACK = STACK(1:end-1);
    XXX = XXX(1:end-1);
end

% find first STACK entry with fit_fields
ii=1;
[~, mod_idxs] = find_modules(STACK, 'nonlinearity', false);
if ~isempty(mod_idxs),
    ii=mod_idxs{end};
    
    [~,mseidx]=find_modules(STACK,'mean_squared_error');
    
    % Then boost on each file individually
    split_stack(ii,mseidx{1}-1);
end

fit05a;


