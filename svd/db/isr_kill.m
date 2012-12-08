% function isr_kill(queueidx,killtype,user)
%
function isr_kill(queueidx,killtype,user)

global DBUSER
if ~exist('user','var') | isempty(user),
   user=DBUSER;
end

if ~exist('killtype','var'),
   killtype=1;
end

if queueidx==-1,
   
   yn=input('kill ALL active queue entries (y/[n])??? ','s');
   if strcmp(yn(1),'y'),
      sql=['SELECT id FROM tQueue',...
           ' WHERE complete=-1 AND user like "',user,'"'];
      queuedata=mysql(sql);
      
      for ii=1:length(queuedata),
         isr_kill(queuedata(ii).id,killtype);
      end
   end
   return
end

sql=['SELECT * FROM tQueue WHERE id=',num2str(queueidx),...
     ' AND user like "',user,'"'];
queuedata=mysql(sql);

if length(queuedata)==0,
   fprintf('ERROR: queueidx=%d not listed in tQueue.\n',queueidx);
   return
elseif queuedata.complete~=-1,
   fprintf('ERROR: queueidx=%d not active.\n',queueidx);
   return
else
   
   cmd=['ssh svd@seil.umd.edu canceljob ',num2str(queuedata.pid)];
   unix(cmd);
   
   fprintf('Killed queueidx=%d (%d).\n',queueidx,killtype);
   dbevent(5,queueidx);
   if killtype==1,
      sql=['UPDATE tQueue set killnow=0 WHERE id=', ...
           num2str(queueidx)];
      mysql(sql);
      isr_start(queueidx);
   else
      sql=['UPDATE tQueue set complete=2,killnow=0 WHERE id=', ...
           num2str(queueidx)];
      mysql(sql);
   end
   
end

