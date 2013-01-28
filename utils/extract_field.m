function ret = extract_field(c, f)
    function v = safeget(x)
        if isfield(x, f)
            v = getfield(x, f);
        else
            v = nan;
        end
    end
    ret = cellfun(@safeget, c, 'UniformOutput', false);
end