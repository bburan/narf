function s = merge_structs(s, s_new)
% Overwrite values of fields in s with values from s_new, if they exist in
% s_new.

if nargin == 2
    fns = fieldnames(s_new);
    for idx = 1:length(fns);
        s.(fns{idx}) = s_new.(fns{idx});
    end
else
    error('merge_structs() must have exactly 2 arguments.');
end