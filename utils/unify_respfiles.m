function xxx = unify_respfiles(xxxs)

xxx = xxxs{1};

for ii = 2:length(xxxs)
    fns = fieldnames(xxxs{ii}.dat);
    for jj = 1:length(fns);
        f = fns{jj};
        xxx.dat.(f) = xxxs{ii}.dat.(f);
    end
    xxx.training_set = union(xxx.training_set, xxxs{ii}.training_set);
    xxx.test_set = union(xxx.test_set, xxxs{ii}.test_set);
    xxx.filecodes = union(xxx.filecodes, xxxs{ii}.filecodes);
    
end

% I'm not sure where the problem is exactly, but this is exploding!
% Do I need to deep copy? Are references building up or something?