function nonNegInteger(property,value)
%nonNegInteger any nonnegative integer
    
%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:26 $

valid =  isreal(value) && isscalar(value) && (value >= 0) && (value == floor(value));
if(~valid)
    error('globaloptim:nonNegInteger:negativeNum','The field ''%s'' must contain a non negative integer.',property);
end
