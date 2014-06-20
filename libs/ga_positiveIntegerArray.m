function positiveIntegerArray(property,value)
%positiveIntegerArray positive integer array

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:39 $

allValid = true;
for i = 1:numel(value)
    valid =  isreal(value(i)) && value(i) == floor(value(i)) && value(i) > 0;
    allValid = allValid && valid;
end
if(~valid)
    error('globaloptim:positiveIntegerArray:notPosIntegerArray','The field ''%s'' must contain a positive integer.',property);
end
