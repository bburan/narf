function xxx = unify_respfiles(xxxs)

xxx = xxxs{1};

for ii = 2:length(xxxs)
    fns = fieldnames(xxxs{ii}.dat);
    for jj = 1:length(fns);
        f = fns{jj};
        xxx.dat.(f) = xxxs{ii}.dat.(f);
    end
end

% TODO: Also need to unify the training_set, test_set, and filecodes!

% I'm not sure where the problem is exactly, but this is exploding!
% Do I need to deep copy? Are references building up or something?