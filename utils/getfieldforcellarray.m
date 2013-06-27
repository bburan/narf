function vals = getfieldforcellarray(c, f)
% vals = getfieldofcellarray(c, f)
% 
% Given a cell array c of structs and a fieldname f, returns a cell array
% of the values found under fieldname f in each struct. If a value is not
% found, a NaN is returned.
%
% ARGUMENTS:
%    c       A cell array of structs
%    f       A field name you want to extract from each struct
%
% RETURNS:
%    vals    A cell array of the values extract from each struct.
% 
    function v = safeget(x)
        if isfield(x, f)
            v = getfield(x, f);
        else
            v = nan;
        end
    end
    vals = cellfun(@safeget, c, 'UniformOutput', false);
end