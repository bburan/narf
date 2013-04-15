function xxx = unify_respfiles(xxxs)

xxx = xxxs{1};

for ii = 2:length(xxxs)
    fns = fieldnames(xxxs{ii}.dat);
    for jj = 1:length(fns);
        f = fns{jj};
        xxx.dat.(f) = xxxs{ii}.dat.(f);
    end
end

