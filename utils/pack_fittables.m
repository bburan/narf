function w = pack_fittables(unused_stack)
% w = pack_fittables()
%
% Packs into a single vector w all the elements of each field of each
% module that has a corresponding field "fit_fields" containing field
% names. If FIT_SPLIT_INDEXES has integers in it, then only params in that
% particular split indexes will be packed. 
%
% Used by fitters to pack and unpack checked values on the STACK. 
%
% Only works for matrix, vector, and numeric types. 

global META STACK;

if ~isfield(META, 'fit_split_indexes')
    only_these_split_indexes = NaN;
else
    only_these_split_indexes = META.fit_split_indexes;
end

w = [];
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
                w = cat(1, w, reshape(m.(p), numel(m.(p)), 1));
            end
        end
    end
end
