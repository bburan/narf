function positiveInteger(property,value)
%positiveInteger any positive integer

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:38 $

valid =  isreal(value) && isscalar(value) && (value > 0) && (value == floor(value));
if(~valid)
   error('globaloptim:positiveInteger:notPosInteger','The field ''%s'' must contain a positive integer.',property);
end

