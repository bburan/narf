% function r=dbsetqueue(queueidx,[progress=0],[complete=-1],[machinename]);
%
% set status for entry in tQueue
% 
% parameters:
% progress - integer value indicating how far along analysis has run
%            (generally while complete=-1)
% complete - 0  = not yet run
%            -1 = in progress
%            1  = complete
%            2  = dead, won't be executed
%
% created SVD 6/2001
%
function r=dbsetqueue(queueidx,progress,complete,machinename);

global USEDB
global DBUSER
global BATQUEUEID

if ~USEDB,
   r=1;
   return;
end

persistent lasttic

if ~exist('queueidx','var') || isempty(queueidx),
   queueidx=BATQUEUEID;
   
   if now-lasttic < 1e-4,
      r=1;
      %now-lasttic
      return
   end
end
lasttic=now;

if isempty(queueidx) | queueidx<=0,
   return
end

dbopen;

if not(exist('complete','var')),
   complete=-1;
end
if not(exist('progress','var')),
   progress=1;
end
if ~exist('machinename','var'),
   
   t=strsep(getenv('MYHOST'),'.');
   if length(t)==0,
      [s,w] = unix('echo $HOSTNAME');
      t=strsep(w,'.');
   end
   
   if s || length(t)==0,
      t{1}='nodeX';
   end
   
   shorthost=deblank(t{1});
   if length(t)>1,
      ext=deblank(t{2});
   else
      ext='isr';
   end
   
   %machinename=[shorthost,'.',ext];
   machinename=shorthost;

   HOSTISLOCAL=1;
else
   HOSTISLOCAL=0;
end

% first grab at queue entry... quickly check to make sure it hasn't
% been grabbed already
if progress==0 & complete==-1,
   sql=['UPDATE tQueue SET',...
        ' complete=',num2str(complete),...
        ' WHERE id=',num2str(queueidx)];
   [res,r]=mysql(sql);
   if r==0,
      fprintf('Queue entry %d just grabbed... Skipping.\n',queueidx);
      return;
   end
end

queuedata=dbgetqueue(queueidx);
if length(queuedata)==0,
   r=0;
   disp('ERROR: dbsetqueue.m:  queue entry not found!');
   if ~isempty(BATQUEUEID) & BATQUEUEID==queueidx,
      BATQUEUEID=[];
   end
   return
end

if (queuedata.pid==0 | progress==0) & complete==-1
   
   % starting job... have to fill in extra fields
   pid=getpid;
   
   computerid=dbgetcomp(machinename);
   if computerid==-1 | isempty(computerid),
      computerid=queuedata.computerid;
   end
   if isempty(computerid),
      computerid=0;
   end
   
   y=year(now);
   mon=month(now);
   d=day(now);
   h=hour(now);
   min=minute(now);
   s=round(second(now));
   startdate=sprintf('%d-%.2d-%.2d %.2d:%.2d:%.2d',y,mon,d,h,min,s);
   
   sql=['UPDATE tQueue SET',...
        ' pid=',num2str(pid),',',...
        ' machinename="',machinename,'",',...
        ' computerid=',num2str(computerid),',',...
        ' startdate="',startdate,'",',...
        ' progress=',num2str(progress),',',...
        ' complete=',num2str(complete),...
        ' WHERE id=',num2str(queueidx)]
   
   [res,r]=mysql(sql);
   
   dbupdateload;
   
   fprintf('Started queue entry %d.\n',queueidx);
   return
end

if ~strcmp(queuedata.user,DBUSER),
   fprintf('Warning: DBUSER=%s setting status for tQueue.user=%s\n',...
           DBUSER,queuedata.user);
end

% check to see if kill flag is set and exit --- unless it's a job
% being managed by dbqueuemaster, in which case the master queuer
% takes care of killing the matlab process
if queuedata.killnow & HOSTISLOCAL,
   fprintf('Queueidx=%d killed (%d).  Exiting Matlab.\n',...
           queueidx,queuedata.killnow);
   if queuedata.killnow==1,
      cstr='0';
   else
      cstr='2';
   end
   sql=['UPDATE tQueue SET',...
        ' killnow=0,',...
        ' complete=',cstr,...
        ' WHERE id=',num2str(queueidx)];
   [res,r]=mysql(sql);
   dbupdateload;
   if queuedata.killnow==3,
      dbdeletequeue(queueidx);
   end
   
   quit
end

% job already exists, just updating entry
%fprintf('Updating job %d (%d).\n',queueidx,progress);
if progress==1,
   progress=queuedata.progress+1;
end

if ~strcmp(queuedata.machinename,machinename) & complete==-1,
   disp('machinename in db doesn''t match current machine! updating db');
   computerid=dbgetcomp(machinename);
   
   sql=['UPDATE tQueue ',...
        'SET progress=',num2str(progress),',',...
        ' machinename="',machinename,'",',...
        ' pid=',num2str(getpid),',',...
        ' computerid=',num2str(computerid),',',...
        ' complete=',num2str(complete),...
        ' WHERE id=',num2str(queueidx)];
else
   sql=['UPDATE tQueue ',...
        'SET progress=',num2str(progress),',',...
        ' complete=',num2str(complete),...
        ' WHERE id=',num2str(queueidx)];
end

[res,r]=mysql(sql);

% always update load to keep dbqueuemaster happy
dbupdateload;

if complete==1,
   dbevent(3);
   
   if ~isempty(queuedata.mailto),
      stitle=sprintf('%s queueid=%d runidx=%d',...
                     queuedata.mailcommand,queuedata.id,queuedata.rundataid);
      mailcommand(queuedata.mailto,...
                  ['disp(''',queuedata.parmstring,'''); ',...
                   queuedata.mailcommand,';'],...
                  stitle);
   end
end


   
