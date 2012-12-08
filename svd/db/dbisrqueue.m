% dbqueuemaster.m:  
%
% control jobs in tQueue, using information in tComputer. runs
% tQueue entries with allowqueuemaster=1 on computers with
% tComputer.allowqueuemaster=1.
% 
% loops until ctrl-break or matlab job is killed. basic flow:
% 1. check for jobs that have tQueue.killnow>0, kill the process on
%    the node and reset job status
% 2. check to see if any jobs or computers have died. if a job has
%    died, mark it as dead to free up space on its computer. if the
%    computer is not responding to pings, remove it from the queue
%    (by setting tComputer.allowqueuemaster=0)
% 3. count number of jobs on each computer
% 4. if a computer is not running at full capacity find the most
%    deserving entry in the queue that hasn't been started (based
%    on number of jobs being run be each user and job priorities)
%    and start it by executing tQueue.progname on the node.
%
% things to be improved:
%
% 1. priorities - betting load balancing. right now things are
% entirely demoncratic, x/n jobs to each of n users for x processors,
% regardless of priority.  all priority does is bump a job to the top
% of the queue for an individual user.  
%    options: kill lower priority jobs with low progress?
%
% 2. overloads/using desktops - if the load on a machine gets too
% high, kill the job running on the machine.  when the load goes
% back down try starting a job again?  would take advantage of
% tComputer.maxload which is currently not used.
%
% 3. dealing with dead computers: it seems that rsh can hang
% sometimes. this makes dbqueuemaster stall without signalling a
% problem.  it has to be killed and restarted to get the queue back
% up
%
% 4. enable facility for sending email when certain jobs have completed.
%
% created SVD 6/9/03
% 6/24/04 SVD - switched from rsh to ssh connection
%

% unix command for running on remote machine via ssh:
SSHCMD=['ssh -q -x -o StrictHostKeyChecking=no',...
        ' -o PasswordAuthentication=no '];

% directory containing shell scripts executed by dbqueuemaster.m
BINPATH='/home/svd/bin/';

% directory containing log files
LOGPATH='/home/tmp/queue/';

% define when "night" starts and stops
SHARESTARTHOUR=21;
SHARESTOPHOUR=8;
NTLOADADJ=0.2;

% establish connection with mysql database
dbopen;

dblog('');
dblog('*************************');
dblog('* DBQUEUEMASTER STARTED *');
dblog('*************************');
dblog('');
dblog('HOST=%s',getenv('MYHOST'));
dblog('BINPATH=%s',BINPATH);
dblog('');

lastreporthour=-1;
lastfairkillmin=minute(now)-10;

%try

% the only way to leave the loop is ctrl-c or killing the matlab process
while 1,
   
   
   %
   % check to see if any jobs have been marked for death
   %
   sql=['SELECT * FROM tQueue WHERE killnow>0 AND allowqueuemaster=1'];
   queuedata=mysql(sql);
   if length(queuedata)>0,
      % kill first job in list
      isr_kill(queuedata(1).id,queuedata(1).killnow);
   end
   
   
   %
   % check to see if jobs or nodes have died
   %
   sql=['SELECT * FROM tQueue WHERE tQueue.complete<0'];
   queuedata=mysql(sql);
   for ii=1:length(queuedata),
      tq=dbgetqueue(queuedata(ii).id);
      if length(tq)>0 & tq.complete==-1 & tq.secago>90+rand*20,
         
         fprintf('check for job %d health\n',tq.id);
         %[s,w]=unix([SSHCMD,'svd@',tq.machinename,...
         %            ' ps h -p ',num2str(tq.pid)]);
         
         if 0& length(w)<3,
            % empty w means that no process with id pid exists
            % record as dead in db
            
            dblog('qid=%d: Proc %d is gone from %s! Setting complete=2.',...
                  tq.id,tq.pid,tq.machinename);
            dbsetqueue(tq.id,tq.progress,2,tq.machinename);
            dbevent(4,tq.id,tq.machinename,tq.pid);
            
            % e-mail notice that job died
            % (if email address is set)
            sql=['SELECT * FROM gUserPrefs WHERE userid="',tq.user,'"'];
            userdata=mysql(sql);
            
            if ~isempty(userdata) & ~isempty(userdata(1).email),
               cmd=['[s,w]=unix(''tail -20 ',LOGPATH,num2str(tq.id),'.out'')'];
               sub=sprintf('Dead job: qid %d (%s) on %s',...
                           tq.id,tq.note,tq.machinename);
               emailres(userdata(1).email,cmd,sub);
            end
         else
            dbupdateload(tq.machinename,tq.user);
            fprintf('queueid=%d: Process %d on %s still alive!\n',...
                    queuedata(ii).id,queuedata(ii).pid,tq.machinename);
            dbset('tQueue',tq.id,'progress',tq.progress+1);
         end
      end
   end
   
   
   % periodically check for dooty emails
   if hour(now)~=lastreporthour | lastreporthour==-1,
      
      calcheck;
      
      lastreporthour=hour(now);
   end
   
   % let things percolate
   pause(1);
   
   unix(''); % flush buffer
end


%catch
   
   % if queue crashes, email someone
   %cmd='fprintf(''Queue crashed!!!\n\n'');'
   %sub='Queue daemon crashed. Need to restart as user queued on nutmeg.';
   %emailres('prenger@socrates.berkeley.edu',cmd,sub);
   

%end
