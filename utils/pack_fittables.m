function w = pack_fittables(~)
% w = pack_fittables()
%
% Packs into a single vector w all the elements of each field of each
% module that has a corresponding field "fit_fields" containing field
% names. 
%
% Used by fitters to pack and unpack checked values on the STACK. 
%
% Only works for matrix, vector, and numeric types. 

global STACK;

w = [];
for ii = 1:length(STACK)
    mm = STACK{ii};
    nsplits = length(mm);    
    for kk = 1:nsplits
        m = mm{kk};        
        if isfield(m, 'fit_fields')
            for jj = 1:length(m.fit_fields),
                p = m.fit_fields{jj};
                w = cat(1, w, reshape(m.(p), numel(m.(p)), 1));
            end
        end
    end
end
