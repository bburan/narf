function names = shorten_modelnames(names)

mnames = {};

for ii = 1:length(names)
    tmp = tokenize_string(names{ii});
    mnames{ii} = [tmp{:}];
end   

% Remove any common prefixes
while length(mnames{1}) > 1 && (all(cellfun(@(a) strcmp(a{1}, mnames{1}{1}), mnames)))
    for ii = 1:length(mnames)
        mnames{ii} = mnames{ii}(2:end);
    end
end

% Remove any common suffixes
while length(mnames{1}) > 1 && (all(cellfun(@(a) strcmp(a{end}, mnames{1}{end}), mnames)))
    for ii = 1:length(mnames)
        mnames{ii} = mnames{ii}(1:end-1);
    end
end

% Flatten the names again
for ii = 1:length(names)
    mm = cellfun(@(x) [x '+'], mnames{ii}, 'UniformOutput', false);
    names{ii} = [mm{:}];
    names{ii} = names{ii}(1:end-1); % Trim off last +, if any    
end
    
end