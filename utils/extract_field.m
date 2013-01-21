function ret = extract_field(c, f)
    ret = cellfun(@(x) getfield(x, f), c, 'UniformOutput', false);
end