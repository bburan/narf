function realUnitScalar(property,value)
%realUnitScalar A scalar on the interval [0,1]

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:45 $

valid = isreal(value) && isscalar(value) && (value >= 0) && (value <= 1);
if(~valid)
    error('globaloptim:realUnitScalar:notScalarOnUnitInterval','The field ''%s'' must contain a scalar on the interval (0,1)',property);
end
