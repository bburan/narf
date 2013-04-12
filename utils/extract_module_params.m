function paramarray = extract_module_params(module_index, paramname) 
% paramarray = extract_module_params(module_index, paramname) 
%
% Returns a cell array of parameters found for the module listed at index
% MODULE_INDEX in the STACK data structure. This function works just fine
% for multiple parameter sets which may exist under STACK. If there is only
% one parameter set, then the returned cell array will have only one
% parameter in it. 
%
% It is intended that you ALWAYS use this function to access data in other
% modules, since you can never know whether there are one or multiple sets
% of parameters for a given module.
%

global STACK;

mm = STACK{module_index};

if iscell(mm)
    paramarray = cell(1, length(mm));
    
    for ii = 1:length(mm)
        m = mm{ii};
        if isfield(m, paramname)
            paramarray{ii} = m.(paramname);
        end
    end
else
    m = mm;
    paramarray = {};
    if isfield(m, paramname)
        paramarray{1} = m.(paramname);
    end
end