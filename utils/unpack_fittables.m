function unpack_fittables(w)
global STACK;

jj = 1;
for ii = 1:length(STACK)
    if isfield(STACK{ii}, 'fit_fields')
        for p = STACK{ii}.fit_fields', p=p{1};
            n = numel(STACK{ii}.(p));
            tmp = w(jj:jj+n-1);
            STACK{ii}.(p) = reshape(tmp, size(STACK{ii}.(p)));
            jj = jj + n;
        end
    end
end