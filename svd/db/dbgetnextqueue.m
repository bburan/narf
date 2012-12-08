% function [BATQUEUEID,queuecomplete,queuedata]=dbgetnextqueue(progname,curbatqueueid,user)
%
% BATQUEUEID -1 means none available
%
% created SVD 4/1/02
%
function [BATQUEUEID,queuecomplete,queuedata]=dbgetnextqueue(progname,...
                                                  curbatqueueid,user);

dbcleanqueue;

global DBUSER
if ~exist('user','var') | isempty(user),
   user=DBUSER;
end

% check to see if max processes are already running on this computer.
[l1,computeridx]=dbupdateload;
sql=['SELECT * FROM tComputer',...
     ' WHERE id=',num2str(computeridx)];
compdata=mysql(sql);
if compdata.numproc>=compdata.maxproc,
   disp('Already at max processes.');
   BATQUEUEID=-1;
   return
end

if ~exist('curbatqueueid','var'),
   BATQUEUEID=0;
else
   BATQUEUEID=curbatqueueid;
end
queuecomplete=-1;
queuedata=[];

while queuecomplete < 0,
   
   sql=['SELECT id,complete FROM tQueue',...
        ' WHERE progname like "',progname,'"',...
        ' AND complete<=0',...
        ' AND user like "',user,'"',...
        ' AND id > ',num2str(BATQUEUEID),...
        ' ORDER BY id '];  % add 'desc' to end of sql to run queue backwards
   batchrundata=mysql(sql);
   
   % find first entry that is not being processed
   rr=1;
   while rr<=length(batchrundata) & batchrundata(rr).complete~=0,
      fprintf('Skipping queue id %d\n',batchrundata(rr).id);
      rr=rr+1;
   end
   
   if rr>length(batchrundata), 
      % no more entries in queue, so quit
      %disp('Queue is empty.');
      BATQUEUEID=-1;
      return
   end
   
   BATQUEUEID=batchrundata(rr).id;
   
   % next check to see if this entry is ready.  someday lock this one
   % as being used at the same time to avoid problems of simultaneous
   % decorr by two processes.
   queuecomplete=dbget('tQueue','complete',BATQUEUEID);
   
   if queuecomplete==0,
      r=dbsetqueue(BATQUEUEID,0,-1);
      if r==0,
         queuecomplete=-1;
      end
   end
   
   if queuecomplete==0,
      % ok, we've confirmed that it's ok, do nothing
      
      sql=['SELECT * FROM tQueue',...
        ' WHERE progname like "',progname,'"',...
        ' AND complete<=0',...
        ' AND user like "',user,'"',...
        ' AND id = ',num2str(BATQUEUEID),...
        ' ORDER BY id '];  % add 'desc' to end of sql to run queue backwards
      queuedata=mysql(sql);
      
      return
      
   elseif queuecomplete==1,
      disp(sprintf('Queue %d decorr already completed.', ...
                   BATQUEUEID));
   elseif queuecomplete==-1,
      disp(sprintf('Queue %d currently running.', ...
                   BATQUEUEID));
   elseif queuecomplete==-2,
      disp(sprintf('Queue %d RC not complete.', ...
                   BATQUEUEID));
   elseif queuecomplete==-3,
      disp(sprintf('Queue %d skipping mail-enabled entry.', ...
                   BATQUEUEID));
   else
      disp(sprintf('Unknown complete code for Queue %d.', ...
                   BATQUEUEID));
   end
end




