% function z = zload(fname,remote)
function z = zload(fname,remote)

if ~exist('remote'),
   remote=0;
end

t = tempname;

if remote,
   if remote==1,
      s=unix(['jcp ',fname,' ',t,'.gz']);
   elseif remote==2,
      s=unix(['bcp ',fname,' ',t,'.gz']);
   end
   
   if s,
      fprintf('can''t find remote file %s\n', fname);
      z=[];
      return
   end
   [s, w] = unix(['gunzip ',t,'.gz']);
else
   
   [s, w] = unix(sprintf('gunzip <%s >%s', fname, t));
   if s,
      fprintf('can''t find %s\n', fname);
      z=[];
      return
   end
end

z=load(t, '-mat');
delete(t);
