function [mdls, xins, xouts] = do_calc_paramsets(stack, xxx)

if ~iscell(stack{end})
    mdls = {stack{end}};
    xins = {xxx{end-1}};
    xouts = {xxx{end}};      
else
    error('Not implemented yet. Should look like core of calc_xxx');
end