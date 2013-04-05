function combos = keyword_combos(mm)
% combos = keyword_combos(mm)
% 
% Returns a cell array of cell arrays of keywords suitable for use with
% fit_single_model(). 
%
% ARGUMENTS:
%   MM       A nested cell array, up to depth level 8. Example:
%            {'env100', 'log2b', {'firn', 'depn'}, {{'npnl', 'mse3'},
%                                                   {'npfnl', 'mse5'}}}
%
% RETURNS:
%   COMBO    A cell array of cell arrays of keywords. Example:
%            {{'env100', 'log2b', 'firn', 'npnl', 'mse3'}, ...
%             {'env100', 'log2b', 'firn', 'npfnl', 'mse5'}, ...
%             {'env100', 'log2b', 'depn', 'npnl', 'mse3'}, ...
%             {'env100', 'log2b', 'depn', 'npfnl', 'mse5'}}

% Single element case
if ischar(mm)
    combos = mm;
    return
end

% Try writing it as tree traversal
% If it's a single element, return it
% If it's a cell array, 
%     prev = {};
%
%     For each element e
%          v = recurse (prev, e)
%          append e to v
% 
%     vs should now be comibined?
%           
%     recurse on each element of the list, providing an index


if iscell(mm) 
    len = zeros(size(mm));
    val = cell(1, length(mm));
    
    for ii = 1:length(mm)
        e = mm{ii};
        ret = keyword_combos(e);
        if ~iscell(ret)
            val{ii} = ret;
        else
            val{ii} = ret{:};
        end 
        
        len(ii) = length(val{ii});
    end
    
    fprintf('mm=%s => ', write_readably(mm));
    fprintf('val=%s\n', write_readably(val));
    
    n = prod(len);
    
    combos = cell(1, n);
    
    for kk = 1:n
        combos{kk} = {};
        inds = [];
        [inds(1), inds(2), inds(3), inds(4), ...
         inds(5), inds(6), inds(7), inds(8)] = ind2sub(len, kk);
        
        for ii = 1:length(mm)
            vi = inds(ii);
            vs = val{ii};
                        
            if ~iscell(vs)
                v = vs;
            else                
                v = vs(vi);
            end
            
            if ~iscell(v)
                combos{kk}{end+1} = v;
            else
                combos{kk} = cat(2, combos{kk}, v{:});
            end
        end
    end
    
    return
end

error('Argument mm was neither a string nor a cell array!');