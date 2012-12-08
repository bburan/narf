function [b,p] = basename(fname)
%
%    func: basename -- return basename of file (ie strip path)
%   usage: basename(fname)
% returns: b == string
%          p == path component
%
%   notes: Mon Aug 11 21:32:15 1997 mazer 
%

v = find(fname == filesep);
if (isempty(v))
  b = fname;
  p = '';
else
  b = fname(v(length(v))+1:length(fname));
  p = fname(1:v(length(v)));
end

