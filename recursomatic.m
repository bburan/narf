function v = recursomatic(o, fn, v)
% o should be a struct with three fields
%    x{1}.depth = 1;
%    x{1}.name = 'A';
%    x{1}.vals = {1 2 3};
    if nargin < 3 v = []; end
    if nargin < 2 error('Recursomatic needs >=2 arguments!'); end
    if isempty(o)
        % TODO: Push v into STACK
        % TODO: Store result somewhere interesting
        RESULT = fn(v);
        return;
    end
    for ii = 1:length(o{1}.vals)
        v.(o{1}.name) = o{1}.vals{ii}; 
        recursomatic(o(2:end), fn, v);
    end
end