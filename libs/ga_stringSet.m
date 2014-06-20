function stringSet(property,value,set)
%stringSet one of a set of strings

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:54 $

for i = 1:length(set)
    if strcmpi(value,set{i})
        return;
    end
end
msg = sprintf('The field %s must contain one of these strings: %s %s %s %s %s',property,set{:});
error('globaloptim:stringSet:notCorrectChoice',msg);
