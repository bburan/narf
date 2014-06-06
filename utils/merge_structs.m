function s = merge_structs(s, s_new)
% s = merge_structs(s, s_new)
%
% Returns a copy of s with all extra fields in s_new included as well. 
% If there is any conflict (ie, the field exists in both), then the value
% from s_new takes precedence. 
%
% ARGUMENTS:
%    s      A structure
%    s_new  A structure
%
% RETURNS:  A modified copy of s

if nargin == 2
    fns = fieldnames(s_new);
    for idx = 1:length(fns);
        if isfield(s, fns{idx}) && isfield(s_new, fns{idx}) && ...
            length(s.(fns{idx})) == 1 && length(s_new.(fns{idx})) == 1 && ...
            isstruct(s.(fns{idx})) && isstruct(s_new.(fns{idx}))
                s.(fns{idx}) = merge_structs(s.(fns{idx}), s_new.(fns{idx}));
        else
            s.(fns{idx}) = s_new.(fns{idx});
        end
    end
else
    error('merge_structs() must have exactly 2 arguments.');
end