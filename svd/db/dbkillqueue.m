% function r=dbkillqueue(queueidx,killtype[=1])
% 
% queueid is the id of an active job (complete=-1) in tQueue
% queueid of -1 kills all jobs in queue (matching current user)
%
% killtype 1 - reset complete to zeros
% killtype 2 - set complete to 2 (don't re-run);
% killtype 3 - kill and delete;
%
function r=dbkillqueue(queueidx,killtype,user)

dbopen;

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
      
      sql=['UPDATE tQueue set killnow=',num2str(killtype),...
           ' WHERE complete=-1',...
           ' AND user like "',user,'"'];
      mysql(sql);
      
      dbqueue(-1);
      
      for ii=1:length(queuedata),
         dbevent(5,queuedata(ii).id);
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
   sql=['UPDATE tQueue set killnow=',num2str(killtype),...
        ' where id=',num2str(queueidx)];
   mysql(sql);
   fprintf('Killed queueidx=%d (%d).\n',queueidx,killtype);
   dbevent(5,queueidx);
end


