function w = pack_fittables(stack)

w = [];
for ii = 1:length(stack)
    if isfield(stack{ii}, 'fit_fields')
        for p = stack{ii}.fit_fields', p=p{1};
            w = cat(1, w, reshape(stack{ii}.(p), numel(stack{ii}.(p)), 1));
        end
    end
end
