function w = pack_fittables(stack)

w = [];
for ii = 1:length(stack)
    if isfield(stack{ii}, 'fit_fields')
        for jj = 1:length(stack{ii}.fit_fields),
            p=stack{ii}.fit_fields{jj};
            w = cat(1, w, reshape(stack{ii}.(p), numel(stack{ii}.(p)), 1));
        end
    end
end
