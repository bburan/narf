function [mdls, xins, xouts] = do_calc_paramsets(stack, xxx)
% do_calc_paramsets 
%
% A helper method

mdls = stack{end};

% TODO: This should look like a splitter!!!!!!
% TODO: This should look more like recalc_xxx()
xins = {xxx{end-1}};
xouts = {xxx{end}};
