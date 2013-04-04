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

kk = 1;
for ii = 1:length(STACK)
    if isfield(STACK{ii}, 'fit_fields')
        for jj = 1:length(STACK{ii}.fit_fields),
            p = STACK{ii}.fit_fields{jj};
            n = numel(STACK{ii}.(p));
            tmp = w(kk:kk+n-1);
            STACK{ii}.(p) = reshape(tmp, size(STACK{ii}.(p)));
            kk = kk + n;
        end
    end
end