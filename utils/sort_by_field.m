function ret = sort_by_field(cellarray, field)
    ff = extract_field(cellarray, field);
    [vals, idxs] = sort(cell2mat(ff));
    ret = cellarray(idxs);
end