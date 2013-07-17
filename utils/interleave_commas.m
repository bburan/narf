function s = interleave_commas(c)
    s = ['"' c{1} '"'];
    for ii = 2:length(c)
        s = cat(2, s, [', "' c{ii} '"']);
    end
end