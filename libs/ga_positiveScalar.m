function positiveScalar(property,value)
%positiveScalar any positive scalar

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:40 $

valid =  isreal(value) && isscalar(value) && (value > 0);
if(~valid)
    error('globaloptim:positiveScalar:notPosScalar','The field ''%s'' must contain a positive scalar.',property);
end
