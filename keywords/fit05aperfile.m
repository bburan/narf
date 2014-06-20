% function fit05aperfile()
%
% use SEMSE rather than standard MSE
function fit05aperfile()

global STACK XXX

% Then boost on each file individually
disp('NOW FITTING EACH FILECODE SEPARATELY');
fit_split_simply(@fit05a, @split_by_filecode, @unify_respfiles); 

