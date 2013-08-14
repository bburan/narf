function unpack_fittables(w)
% unpack_fittables(w, only_these_split_indexes)
%
% Given a weight wector w, this function fills up the current fields marked
% in the STACK as being fittable with the values from w. It is the opposite
% action to pack_fittables(). If ONLY_THIS_SPLIT_INDEX is an integer, then only values in that
% particular split indexes will be unpacked. 
%
% ARGUMENTS:
%    w       A weight vector as created by pack_fittables().
%
% RETURNS: Nothing

global STACK META;

if ~isfield(META, 'fit_split_indexes')
    only_these_split_indexes = NaN;
else
    only_these_split_indexes = META.fit_split_indexes;
end

rr = 1;
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
                n = numel(m.(p));
                tmp = w(rr:rr+n-1);
                STACK{ii}{kk}.(p) = reshape(tmp, size(m.(p)));
                rr = rr + n;
            end
        end
    end
end