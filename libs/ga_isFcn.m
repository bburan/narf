function [handle,args] =  isFcn(fun)
%isFcn    get function handle and args (if fun is in cell array form)    

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:26:22 $

% If it's a scalar fcn handle or a cellarray starting with a fcn handle and
% followed by something other than a fcn handle, return parts, else empty    
  handle = [];
  % If fun is a cell array with additional arguments, handle them
  if iscell(fun) 
      if ~isempty(fun)
          args = fun(2:end);
          handle = fun{1};
      else  % Cell could be empty too
          args = {};
      end
  else % Not a cell
      args = {};
      handle = fun;
  end
  
  if ~isempty(handle)
      [handle,msg] = fcnchk(handle);
      if ~isempty(msg)
          handle =[];
      end
  end