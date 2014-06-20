function [handle,args] = ga_functionHandleOrCell(property,value)
%functionHandleOrCell A function Handle or a cell array starting with a function
%handle.

%   Copyright 2007-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:10 $

[handle,args] = ga_isFcn(value);

if ~isempty(handle)
    return
elseif strcmp(property,'NonconFcn')
    error('globaloptim:functionHandleOrCell:needFunctionHandle','The constraint function must be a function handle.');
else
    error('globaloptim:functionHandleOrCell:needFunctionHandle','The field ''%s'' must contain a function handle.',property);
end

