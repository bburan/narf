function [mdls, xins, xouts] = calc_paramsets(stack, xxx)
% [mdls, xins, xouts] = calc_paramsets(stack, xxx)
%
% Used by many module plot routines. Just makes a convenient form of
% copies of data. 
%
% RETURNS
%    mdls   Actually just stack{end}, an array of module paramsets
%    xins   Cell array of inputs to each of the module paramsets
%    xouts  Cell array of the outputs from each of the module paramsets

mdls = stack{end};

calc_xxx(length(stack), length(stack));

xins = cell(1, length(mdls));
xouts = cell(1, length(mdls));

for jj = 1:length(mdls)
    xins{jj} = xxx{end-1};
    xouts{jj} = xxx{end};
end
