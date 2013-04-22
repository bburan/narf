function [mdls, xins, xouts] = calc_paramsets(stack, xxx)
% [mdls, xins, xouts] = calc_paramsets(stack, xxx)
%
% Used by calc_xxx() and many module plot routines to iterate through each
% parameter set in the end of the STACK (i.e. STACK{end}) and compute the
% inputs and outputs to each, before the unifier function collapses the
% data into a single. Therefore, this function is useful when you need to 
% get at the data that would otherwise be hidden by the unifier's 
% collapsing action. 
%
% RETURNS
%    mdls   Actually just stack{end}, an array of module paramsets
%    xins   Cell array of inputs to each of the module paramsets
%    xouts  Cell array of the outputs from each of the module paramsets

mdls = stack{end};

if isfield(mdls{1}, 'splitter')
    xins = mdls{1}.splitter(xxx);
else
    xins = cell(size(mdls{1}));
    xins(:) = {xxx};
end

if length(xins) ~= length(mdls)
    error('Splitter group count does not match STACK parameter set count.');
end

if ~mdls{1}.isready_pred(stack, xxx);
    error('Stack was not fully ready for module %s at index %d', ...
          mdls{1}.name, length(stack));
end

xouts = cell(1, length(mdls));

% Iterate through each parameter set if there is a unifier
for jj = 1:length(mdls)
    xouts{jj} = mdls{jj}.fn(mdls{jj}, xins{jj}{end}, stack, xins);
end
