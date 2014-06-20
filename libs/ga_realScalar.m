function realScalar(property,value)
%realScalar Test for real scalar

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:44 $

valid = isreal(value) && isscalar(value);
if(~valid)
    error('globaloptim:realScalar:notScalar','The field ''%s'' must contain a scalar',property);
end
