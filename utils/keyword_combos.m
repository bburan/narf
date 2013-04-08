function combos = keyword_combos(mm)
% combos = keyword_combos(mm)
% 
% Returns a cell array of cell arrays of keywords suitable for use with
% fit_single_model(). 
%
% ARGUMENTS:
%   MM       A nested cell array. Example:
%            {'env100', 'log2b', {'firn', 'depn'}, {{'npnl', 'mse3'},
%                                                   {'npfnl', 'mse5'}}}
%
% RETURNS:
%   COMBO    A cell array of cell arrays of keywords. Example:
%            {{'env100', 'log2b', 'firn', 'npnl', 'mse3'}, ...
%             {'env100', 'log2b', 'firn', 'npfnl', 'mse5'}, ...
%             {'env100', 'log2b', 'depn', 'npnl', 'mse3'}, ...
%             {'env100', 'log2b', 'depn', 'npfnl', 'mse5'}}

paths = {};

function traverse (tree, path)
    %fprintf('traverse(%s,\n\t %s)\n\n', write_readably(tree), write_readably(path));
    
    if isempty(tree)    % On a leaf
        paths{end+1} = path;   
        
    elseif iscell(tree) % On an intermediate node        
        head = tree{1};
        tail = tree(2:end);
        if ischar(head)
            traverse(tail, cat(2, path, head));
        elseif iscell(head)
            for ii = 1:length(head)
                traverse(tail, cat(2, path, head{ii}));
            end
        else
            error('Neither a string nor a cell array: %s', write_readably(tree));
        end
    else
        error('Neither a string nor a cell array: %s', write_readably(tree));    
    end
end

traverse(mm, {});

combos = paths;

end