function unpack_fittables(w)
% unpack_fittables(w)
%
% Given a weight wector w, this function fills up the current fields marked
% in the STACK as being fittable with the values from w. It is the opposite
% action to pack_fittables().
%
% ARGUMENTS:
%    w       A weight vector as created by pack_fittables().
%
% RETURNS: Nothing

global STACK;

rr = 1;
for ii = 1:length(STACK)    
    mm = STACK{ii};
    
    for kk = 1:length(mm)
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