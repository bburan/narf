function ret = only_named_summaries(summaries, modelnames)
% Returns a cell array of elements of SUMMARIES whose modelnames are also
% found in cell array MODELNAMES. 

ret = {};
for ii = 1:length(summaries)
    s = summaries{ii};
    if any(cellfun(@(m)isequal(s.modelname, m), modelnames))
        ret{end+1} = s;
    end
end



