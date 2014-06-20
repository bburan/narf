function nonNegScalar(property,value)
%nonNegScalar any scalar >= 0    

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:27 $

valid =  isreal(value) && isscalar(value) && (value >= 0);
if(~valid)
    error('globaloptim:nonNegScalar:notNonNegScalar','The field ''%s'' must contain a non negative scalar.',property);
end
