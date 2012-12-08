% isr_queuemaster.m:  
%
% much-reduced functionality from dbqueuemaster.m, since the isr
% cluster does most of the work now.
%
% main jobs:
% 1. check to see if quiet jobs are still alive. mark dead if not
% 2. check to see if any jobs are marked for death. if so, tell the
%    cluster to kill them
% 3. check to see if anyone needs to be emailed.
%
% created SVD 2007-05-30
%

% unix command for running on remote machine via ssh:
SSHCMD=['ssh -q -x -o StrictHostKeyChecking=no',...
        ' -o PasswordAuthentication=no '];

% directory containing shell scripts executed by dbqueuemaster.m
BINPATH='/auto/p1/svd/bin/';

% directory containing log files
LOGPATH='/auto/data/tmp/queue/';

% define when "night" starts and stops
SHARESTARTHOUR=21;
SHARESTOPHOUR=8;
NTLOADADJ=0.2;

% establish connection with mysql database
dbopen;

dblog('');
dblog('***************************');
dblog('* ISR_QUEUEMASTER STARTED *');
dblog('***************************');
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
   % check to see if there are any jobs to start
   %
   sql=['SELECT * FROM tQueue WHERE complete=0 ORDER BY id'];
   queuedata=mysql(sql);
   if length(queuedata)>0,
      % start first job in list
      isr_start(queuedata(1).id);
      plen=1+rand*3;
      fprintf('pausing %.1f seconds to percolate...\n',plen);
      pause(plen);
   end
   
   
   %
   % check to see if any jobs have been marked for death
   %
   sql=['SELECT * FROM tQueue WHERE killnow>0'];
   queuedata=mysql(sql);
   if length(queuedata)>0,
      % kill first job in list
      isr_kill(queuedata(1).id,queuedata(1).killnow);
   end
   
   
   %
   % check to see if jobs or nodes have died
   %
   sql=['SELECT *,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
        ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
        ' FROM tQueue WHERE tQueue.complete<0'];
   queuedata=mysql(sql);
   for ii=1:length(queuedata),
      tq=queuedata(ii);
      if length(tq)>0 & tq.complete==-1 & tq.secago>240+rand*20 & tq.pid>0,
         
         [status,extended]=isr_checkjob(tq.id);
         
         if status,
            fprintf('queueid=%d: Process %d on %s still alive!\n',...
                    queuedata(ii).id,queuedata(ii).pid,tq.machinename);
            dbset('tQueue',tq.id,'progress',tq.progress+1);
         else
            tq=dbgetqueue(queuedata(ii).id);
            
            if tq.complete==-1,
               
               dblog('qid=%d: Proc %d is gone from %s! Setting complete=2.',...
                     tq.id,tq.pid,tq.machinename);
               dbsetqueue(tq.id,tq.progress,2);
               dbevent(4,tq.id,tq.machinename,tq.pid);
               
               % e-mail notice that job died
               % (if email address is set)
               sql=['SELECT * FROM gUserPrefs WHERE userid="',tq.user,'"'];
               userdata=mysql(sql);
               
               if 0& ~isempty(userdata) & ~isempty(userdata(1).email),
                  cmd=['[s,w]=unix(''tail -20 ',LOGPATH,num2str(tq.id),'.out'')'];
                  sub=sprintf('Dead job: qid %d (%s) on %s',...
                              tq.id,tq.note,tq.machinename);
                  emailres(userdata(1).email,cmd,sub);
               end
            end
         end
      end
   end
   
   
   % periodically check for dooty emails
   if hour(now)~=lastreporthour | lastreporthour==-1,
      
      calcheck;
      
      lastreporthour=hour(now);
   end
   
   % tell the database dbqueuemaster is still alive
   [s,dhost]=unix('hostname');
   dhost=strtrim(dhost);
   sql=['UPDATE tGlobalData SET daemonclick=now(),daemonhost="',dhost,'"'];
   mysql(sql);
   
   % let things percolate
   pause(0.2);
   
   [s,w]=unix(''); % flush buffer
end


%catch
   
   % if queue crashes, email someone
   %cmd='fprintf(''Queue crashed!!!\n\n'');'
   %sub='Queue daemon crashed. Need to restart as user queued on nutmeg.';
   %emailres('prenger@socrates.berkeley.edu',cmd,sub);
   

%end
