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
    mm = stack{ii};
    
    if iscell(mm)        
        for kk = 1:length(mm)
            m = mm{kk};
            
            if isfield(m, 'fit_fields')
                for jj = 1:length(m.fit_fields),
                    p = m.fit_fields{jj};
                    w = cat(1, w, reshape(m.(p), numel(m.(p)), 1));
                end
            end
        end
    else
        m = mm;
        if isfield(stack{ii}, 'fit_fields')
            for jj = 1:length(stack{ii}.fit_fields),
                p = stack{ii}.fit_fields{jj};
                w = cat(1, w, reshape(stack{ii}.(p), numel(m.(p)), 1));
            end
        end
    end
end
