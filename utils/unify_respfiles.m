function xxx = unify_respfiles(xxxs)
 
xxx = xxxs{1};
%keyboard
for ii = 2:length(xxxs)
    fns = fieldnames(xxxs{ii}.dat);
    for jj = 1:length(fns);
        f = fns{jj};
        xxx.dat.(f) = xxxs{ii}.dat.(f);
    end
    % union sorts by alphabetical order and can cause mismatch in order of
    % file names and filecodes:
    %xxx.training_set = union(xxx.training_set, xxxs{ii}.training_set);
    %xxx.test_set = union(xxx.test_set, xxxs{ii}.test_set);
    %xxx.filecodes = union(xxx.filecodes, xxxs{ii}.filecodes);
    xxx.training_set = cat(2,xxx.training_set, xxxs{ii}.training_set);
    xxx.test_set = cat(2,xxx.test_set, xxxs{ii}.test_set);
    xxx.filecodes = cat(2,xxx.filecodes, xxxs{ii}.filecodes);
end

% I'm not sure where the problem is exactly, but this is exploding!
% Do I need to deep copy? Are references building up or something?