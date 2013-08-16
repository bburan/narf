function [depths]= find_fit_param_depths(stack)
% [depths] = find_fit_param_depths(stack)
%
% Returns a vector of indexes of STACK containing free parameters, ordered
% in the same manner as pack_fittables() and unpack_fittables().
%
% This vector can be used to reduce the amount of work that fitting
% algorithms do, since they only need to recalculate the stack from the 
% point where fields have changed since the previous loop. 

global META STACK;

if ~isfield(META, 'fit_split_indexes')
    only_these_split_indexes = NaN;
else
    only_these_split_indexes = META.fit_split_indexes;
end

depths = [];
for ii = 1:length(STACK)
    mm = STACK{ii};
    nsplits = length(mm);    
    for kk = 1:nsplits
        % If this is a split module but not in the right indexes, skip it
        if nsplits > 1 && ~isnan(only_these_split_indexes) && ...
               ~isempty(find(kk == only_these_split_indexes, 1)),  
            continue;
        end
        m = mm{kk};        
        if isfield(m, 'fit_fields')
            for jj = 1:length(m.fit_fields),
                p = m.fit_fields{jj};
                w = reshape(m.(p), numel(m.(p)), 1);
                depths = cat(1, depths, ii*ones(size(w)));
            end
        end
    end
end
