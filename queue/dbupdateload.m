% function [l1,computerid]=dbupdateload(hostname,user,updateloadbg);
%
% read load for hostname from /proc/loadavg and update the info in
% tComputer (for access via dbgetload and dbloadmon)
%
% for any computer on the jlg.berkeley.edu subnet, also checks to
% see if the processes listed in the active tQueue entries are
% actually alive. if not, sets them to complete=0
%
% created SVD 6/01
% modified SVD 1/29/03 - added check to see if jobs are still alive
%
function [l1,computerid]=dbupdateload(hostname,user,updateloadbg);

SSHCMD=['ssh -q -x -o ConnectTimeout=3 -o StrictHostKeyChecking=no',...
        ' -o PasswordAuthentication=no '];

dbopen;
if ~exist('hostname','var'),
   hostname=getenv('MYHOST');
   if isempty(hostname),
      [s,hostname]=unix('hostname');
      if ~isempty(find(hostname==' ' | hostname==char(10))),
         [s,hostname]=unix('hostname');
      end
      hostname=deblank(hostname);
   end
   hostislocal=1;
else
   hostislocal=0;
end

if ~exist('user','var'),
   user=getenv('USER');
end

if strcmp(user,'root') | strcmp(user,''),
   user='svd';
end
if ~exist('updateloadbg','var'),
   updateloadbg=0;
end
if updateloadbg,
   bgstring=' &';
else
   bgstring='';
end

[computerid,compdata]=dbgetcomp(hostname);
if computerid==-1,
   return
end

% check how many active processes there are on the machine,
% according to the database
sql=['SELECT * from tQueue',...
     ' WHERE complete=-1',...
     ' AND computerid=',num2str(computerid)];
queuedata=mysql(sql);
numproc=length(queuedata);

% get load either locally or via ssh
if hostislocal,
   r=load('/proc/loadavg');
   
   l1=r(1);
   l5=r(2);
   l15=r(3);
   
   % update the number of jobs in the database
   if numproc~=compdata.numproc,
      sql=['UPDATE LOW_PRIORITY tComputer SET',...
           ' load1=',num2str(l1),',',...
           ' load5=',num2str(l5),',',...
           ' load15=',num2str(l15),',',...
           ' pingcount=0,',...
           ' numproc=',num2str(numproc),...
           ' WHERE id=',num2str(computerid)];
      [r,aff]=mysql(sql);
   end
else
   sql=['UPDATE LOW_PRIORITY tComputer SET',...
        ' pingcount=pingcount+1,',...
        ' numproc=',num2str(numproc),...
        ' WHERE id=',num2str(computerid)];
   [r,aff]=mysql(sql);
   
   shorthostname=strsep(hostname,'.');
   shorthostname=shorthostname{1};
   tuser='svd@'; % [user '@'];
   %disp([SSHCMD,tuser,shorthostname,' /auto/users/svd/bin/qsetload',bgstring]);
   [s,w]=unix([SSHCMD,tuser,shorthostname,' /auto/users/svd/bin/qsetload',bgstring]);
   
   if s>0,
      fprintf('dbupdateload.m: s=%d re-checking load on: %s\n',s,hostname);
      [s,w]=unix([SSHCMD,tuser,shorthostname,' /auto/users/svd/bin/qsetload']);
      fprintf('failed with s=%d w=%s\n',s,w);
   end
   
   % if failed, return -1. this lets dbqueuemaster know that the
   % computer is not responding.
   if s>0,
      l1=-1;
   else 
      % get info about hostname
      l1=compdata.load1;
   end

   % update the number of jobs in the database
   if numproc~=compdata.numproc,
      sql=['UPDATE LOW_PRIORITY tComputer SET',...
           ' numproc=',num2str(numproc),...
           ' WHERE id=',num2str(computerid)];
      [r,aff]=mysql(sql);
   end
end

return

%skip this last thing because pids are not accurate!


% for jobs not running via dbqueuemaster, check to see if they're
% still alive.
for nn=1:numproc,
   
   % check to see if queue process still exists on hostname
   % if allowqueuemaster, dbqueuemaster takes care of this
   if ~queuedata(nn).allowqueuemaster,
      
      pid=queuedata(nn).pid;
      if hostislocal,
         [s,w]=unix(['ps h -p ',num2str(pid)]);
      elseif strcmp(compdata.ext,'jlg') | strcmp(compdata.ext,'fet'),
         [s,w]=unix(['rsh -l david ',hostname,' ps h -p ',num2str(pid)]);
      elseif strcmp(compdata.ext,'isr'),
         ['ssh -l svd ',hostname,' ps h -p ',num2str(pid)]
         [s,w]=unix(['ssh -l svd ',hostname,' ps h -p ',num2str(pid)]);
      else
         w='xxxx'; % force a skip, can't update remotely.
      end
      if length(w)<3,
         
         % process doesn't exist! record as dead in db
         fprintf('Process %d on %s is dead! Reseting to complete=0.\n',...
                 pid,hostname);
         dbsetqueue(queuedata(nn).id,0,0);
         numproc=numproc-1;
      end
   end 
end
%fprintf('host: %s loc: %s numproc: %d load: %.2f %.2f %.2f\n',...
%        shorthost,ext,numproc,l1,l5,l15);
