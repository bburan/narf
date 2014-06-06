function [depths]= find_fit_param_depths(~)
% [depths] = find_fit_param_depths()
%
% Returns a vector of indexes of STACK containing free parameters, ordered
% in the same manner as pack_fittables() and unpack_fittables().
%
% This vector can be used to reduce the amount of work that fitting
% algorithms do, since they only need to recalculate the stack from the 
% point where fields have changed since the previous loop. 

global STACK;

depths = [];
depth = 1;
for ii = 1:length(STACK)
    mm = STACK{ii};
    nsplits = length(mm);    
    for kk = 1:nsplits
        m = mm{kk};        
        if isfield(m, 'fit_fields')
            for jj = 1:length(m.fit_fields),
                p = m.fit_fields{jj};
                w = reshape(m.(p), numel(m.(p)), 1);
                depths = cat(1, depths, depth*ones(size(w)));
            end
        end
        depth = depth + 1;
    end
end
