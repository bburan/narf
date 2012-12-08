function z = zload(fname, ascii)
%function z = zload(fname, ascii)
%
% Load a (possibly) gziped file into matlab as data.  If ascii is
% true, then treat as ascii table file.  Otherwise, assume this
% is a matrix.  If a return value is given, then the loaded matrix
% or struct is returned.  Otherwise, load occurs in caller's
% stack frame.
%
% Mon Apr  9 09:30:23 2001 mazer 

if ~exist('ascii', 'var')
  ascii = 0;
end

if ascii
  type = '-ascii';
else
  type = '-mat';
end

%fname = ex(fname);
  
if gzipped(fname) & ~exist(fname,'file'),
   t=strrep(fname,'.gz','');
   cleanup = 0;
   
elseif gzipped(fname),
   t = tempname;
   [s, w] = unix(sprintf('gunzip <%s >%s', fname, t));
   if s
      error(sprintf('can''t find %s', fname));
   end
   cleanup = 1;
else
   t = fname;
   cleanup = 0;
end

if nargout == 0
  evalin('caller', sprintf('load(''%s'', ''%s'')', t, type));
else
  z = load(t, type);
end

if cleanup
  delete(t);
end

function f = ex(varargin)
%function p = ex(varargin)
%
%  tilde expand a filename..
%
% Tue Feb  6 14:00:06 2001 mazer 

f = '';
for n=1:length(varargin)
  s = varargin{n};
  f = [f strrep(s, '~', getenv('HOME'))];
end


