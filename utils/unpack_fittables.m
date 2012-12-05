function unpack_fittables(w)
global STACK;

jj = 1;
for ii = 1:length(STACK)
    if isfield(STACK{ii}, 'fit_fields')
        for p = STACK{ii}.fit_fields', p=p{1};
            n = numel(STACK{ii}.(p));
            tmp = w(jj:n);
            STACK{ii}.(p) = reshape(tmp, size(STACK{ii}.(p)));
            j =+ n;
        end
    end
end