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
        s.(fns{idx}) = s_new.(fns{idx});
    end
else
    error('merge_structs() must have exactly 2 arguments.');
end