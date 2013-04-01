function w = pack_fittables(stack)
% w = pack_fittables(stack)
%
% Packs into a single vector w all the elements of each field of each
% module that has a corresponding field "fit_fields" containing field
% names. Got it? ;-)
%
% Used by fitters to pack and unpack checked values on the stack. 
%
% Only works for matrix, vector, and numeric types. 

w = [];
for ii = 1:length(stack)
    if isfield(stack{ii}, 'fit_fields')
        for jj = 1:length(stack{ii}.fit_fields),
            p = stack{ii}.fit_fields{jj};
            w = cat(1, w, reshape(stack{ii}.(p), numel(stack{ii}.(p)), 1));
        end
    end
end
