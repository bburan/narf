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

narf_set_path

% unix command for running on remote machine via ssh:
SSHCMD=['ssh -q -x -o StrictHostKeyChecking=no',...
        ' -o PasswordAuthentication=no '];
SSHCMD=['sudo su %s -c ''ssh -q -x -o ConnectTimeout=3 -o StrictHostKeyChecking=no',...
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

lastsudomin=-1;
lastreporthour=-1;
lastfairkillmin=minute(now)-10;

%try

% the only way to leave the loop is ctrl-c or killing the matlab process
while 1,
   
   % periodically make sure sudo is validated
   if minute(now)~=lastsudomin,
       !sudo ls
       lastsudomin=minute(now);
   end
   
   %
   % check to see if any jobs need to be killed on night-only machines 
   %
   nighttime=(hour(now)>=SHARESTARTHOUR | hour(now)<SHARESTOPHOUR);
   if ~nighttime,
      sql=['SELECT tQueue.* FROM tComputer,tQueue',...
           ' WHERE tComputer.id=tQueue.computerid',...
           ' AND tQueue.complete=-1',...
           ' AND tComputer.dead=0',...
           ' AND tComputer.allowqueuemaster=3'];
      tq=mysql(sql);
      if length(tq)>0,
         dbset('tQueue',tq(1).id,'killnow',1);
         dblog('node %s: night is over. killing qid %d',...
               tq(1).machinename,tq(1).id);
         dbevent(5,tq(1).id,tq(1).machinename,tq(1).pid,tq(1).user);
      end
   end
   
   %
   % check to see if any machines are overloaded
   % (allowqueuemaster=2 and load>killqueueload) and mark their
   % jobs for death
   %
   if nighttime,
      nextra=['+',num2str(NTLOADADJ)];
   else
      nextra='';
   end
   sql=['SELECT tQueue.*,tComputer.killqueueload FROM tComputer,tQueue',...
        ' WHERE tComputer.id=tQueue.computerid',...
        ' AND tComputer.allowqueuemaster in (2,3)',...
        ' AND tComputer.load1>(tComputer.killqueueload',nextra,')',...
        ' AND tQueue.complete=-1'];
   tq=mysql(sql);
   if length(tq)>0,
      % only kill one job. might make load go down enough to keep
      % any others that are running (for machines with maxproc>1)
      dbset('tQueue',tq(1).id,'killnow',1);
      dbset('tComputer',tq(1).computerid,'lastoverload',1);
      dblog('qid %d: killing - node %s over load %.2f. ',...
            tq(1).id,tq(1).machinename,...
            tq(1).killqueueload+nighttime*NTLOADADJ);
      dbevent(5,tq(1).id,tq(1).machinename,tq(1).pid,tq(1).user);
   end
   
   %
   % check to see if any machines marked as overloaded are
   % available, and mark them as not overloaded
   %
   sql=['SELECT * FROM tComputer',...
        ' WHERE tComputer.lastoverload=1',...
        ' AND tComputer.dead=0',...
        ' AND load1<allowqueueload',nextra,...
        ' AND load15<allowqueueload',nextra];
   compdata=mysql(sql);
   for ii=1:length(compdata),
      fprintf('reactivating cond machine %s\n',compdata(ii).name);
      dbset('tComputer',compdata(ii).id,'lastoverload',0);
   end
      
   %sql=['SELECT * ,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
   % ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
   % ' FROM tComputer order by load1,name'];
   
   %
   % check to see if any jobs have been marked for death
   %
   sql=['SELECT * FROM tQueue WHERE killnow>0 AND allowqueuemaster=1'];
   queuedata=mysql(sql);
   if length(queuedata)>0,
      
      % kill first job in list
      cmd=sprintf([SSHCMD,queuedata(1).machinename,...
            ' "',BINPATH,'childrenof ',num2str(queuedata(1).pid),...
            ' | xargs kill"'''],queuedata(1).user);
      %cmd=[SSHCMD,queuedata(1).user,'@',queuedata(1).machinename,...
      %           ' "',BINPATH,'childrenof ',num2str(queuedata(1).pid),...
      %            ' | xargs kill" > /dev/null &'];
      disp(cmd);
      [s,w]=unix(cmd);
      %[s,w]=unix([SSHCMD,queuedata(1).user,'@',queuedata(1).machinename,...
      %           ' "',BINPATH,'childrenof ',num2str(queuedata(1).pid),...
      %            ' | xargs kill" > /dev/null &']);
      
      % clean up temp directory on the machine where the job was killed
      cmd=sprintf([SSHCMD,queuedata(1).machinename,...
                   ' "\rm -R /tmp/',num2str(queuedata(1).id),...
                   '" > /dev/null &'''],...
                  queuedata(1).user);
      %cmd=[SSHCMD,queuedata(1).user,'@',queuedata(1).machinename,...
      %        ' "\rm -R /tmp/',num2str(queuedata(1).id),'" > /dev/null &']
      disp(cmd);
      [s,w]=unix(cmd);
      
      tq=dbgetqueue(queuedata(1).id);
      
      if tq.killnow>0,
         
         dbset('tQueue',tq.id,'killnow',0);
         
         if tq.killnow==2,
            dblog('qid=%d: received kill & leave dead signal (%s pid=%d)',...
                  tq.id,tq.machinename,tq.pid);
            dbsetqueue(tq.id,tq.progress,2,tq.machinename);
         elseif tq.killnow==3,
            dblog('qid=%d: received kill & delete signal (%s pid=%d)',...
                  tq.id,tq.machinename,tq.pid);
            dbsetqueue(tq.id,tq.progress,2,tq.machinename);
            dbdeletequeue(tq.id,tq.user);
         else
            dblog('qid=%d: received kill signal (%s pid=%d)',...
                  tq.id,tq.machinename,tq.pid);
            dbsetqueue(tq.id,0,0,tq.machinename);
         end
      end
   end
   
   %
   % check to see if jobs or nodes have died
   %
   sql=['SELECT * FROM tQueue WHERE tQueue.complete<0'];
   queuedata=mysql(sql);
   for ii=1:length(queuedata),
      tq=dbgetqueue(queuedata(ii).id);
      if length(tq)>0 & tq.complete==-1 & tq.secago>90+rand*20,
         
          CMD=sprintf([SSHCMD,tq.machinename,...
                       ' ps h -p ',num2str(tq.pid),''''],'svd');
          [s,w]=unix(CMD);
          
          % recheck status to make sure that the job didn't JUST
          % complete (which happens, esp on a dbkillqueue(-1)
          tq=dbgetqueue(tq.id);
          
          if s==255,
              
              % ssh failed... mark computer as dead
              dblog('host=%s: ssh check failed, removing from cluster',...
                    tq.machinename);
              dbset('tComputer',tq.computerid,'dead',1);
              dbsetqueue(tq.id,0,0,tq.machinename);
              dbevent(7,tq.id,tq.machinename);
              dbevent(5,tq.id,tq.machinename,tq.pid);
              
          elseif tq.complete==-1 & length(w)<3,
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
                  disp('skipping email');
                  %emailres(userdata(1).email,cmd,sub);
              end
          else
            dbupdateload(tq.machinename,tq.user);
            fprintf('queueid=%d: Process %d on %s still alive!\n',...
                    queuedata(ii).id,queuedata(ii).pid,tq.machinename);
            dbset('tQueue',tq.id,'progress',tq.progress+1);
          end
      end
   end

   
   %
   % make sure load measures and job counts are up-to-date
   %
   
   % first update computer load
   %fprintf('ldchk: <');
   %dbgetload([],1);
   dbgetload([],1,20,1);
   %fprintf('>');
   
   % count active jobs in tQueue and update numproc in tComputer if
   % it's wrong
   sql=['SELECT count(tQueue.id) as qnumproc,',...
        ' tComputer.id,tComputer.name,tComputer.numproc',...
        ' FROM tComputer LEFT JOIN tQueue',...
        ' ON tQueue.computerid=tComputer.id',...
        ' WHERE tQueue.complete=-1',...
        ' GROUP BY tComputer.id',...
        ' HAVING qnumproc<>numproc'];
   %     ' WHERE tComputer.allowqueuemaster>0',...
   jobcountdata=mysql(sql);
   
   for ii=1:length(jobcountdata),
      if jobcountdata(ii).qnumproc~=jobcountdata(ii).numproc
         sql=['UPDATE tComputer',...
              ' SET numproc=',num2str(jobcountdata(ii).qnumproc),...
              ' WHERE id=',num2str(jobcountdata(ii).id)];
         mysql(sql);
      end
   end
   
   sql=['SELECT * FROM tComputer WHERE dead=0 AND pingcount>4'];
   compdata=mysql(sql);
   for ii=1:length(compdata),
      dblog('host=%s.%s: pingcount>4 dead? removing from cluster',...
            compdata(ii).name,compdata(ii).ext);
      dbset('tComputer',compdata(ii).id,'dead',1);
      dbset('tComputer',compdata(ii).id,'pingcount',0);
      dbevent(7,0,compdata(ii).id);
   end
   
   sql=['SELECT * FROM tComputer WHERE dead',...
        ' AND ((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
        ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate))>3600'];
   compdata=mysql(sql);
   for ii=1:length(compdata),
      dblog('host=%s.%s: dead 6 hours; attempting to resurrect',...
            compdata(ii).name,compdata(ii).ext);
      dbset('tComputer',compdata(ii).id,'dead',0);
      dbset('tComputer',compdata(ii).id,'pingcount',0);
   end
   
   %
   % find machines with fewer than maximum jobs
   %
   % query optimized using info from previous queries in this loop
   % of the daemon.
   %
   if nighttime,
      sallow='allowqueuemaster>0';
   else
      sallow='allowqueuemaster in (1,2)';
   end
   sql=['SELECT * FROM tComputer',...
        ' WHERE numproc<maxproc AND ',sallow,...
        ' AND lastoverload=0 AND dead=0',...
        ' ORDER BY allowqueuemaster,load1'];
   compdata=mysql(sql);
   
   if length(compdata)==0,
       %disp('no computers available, pausing');
      pause(1);
   end
   nextid=0;  % flag 0 to avoid "not at capacity" message
   
   compidx=0;
   compok=0;
   while compidx<length(compdata) & ~compok,
      compidx=compidx+1;
      
      % check to see whether this computer is ok for giving a job
      c=compdata(compidx);
      if c.allowqueuemaster==1,
         compok=1;
      elseif c.allowqueuemaster==2 & (~c.lastoverload & ...
                       c.load1<c.killqueueload+NTLOADADJ*nighttime),
         % machine has to not be overloaded and below auto-kill load
         compok=1;
      elseif c.allowqueuemaster==3 & nighttime & ...
             (~c.lastoverload & c.load1<(c.killqueueload+NTLOADADJ)),
         % machine allows jobs at night
         fprintf('Yee haw! Night time! Cond maxload +%.1f\n',NTLOADADJ);
         compok=1;
      end
      
      % find the most deserving user for the next computer in the queue
      if c.allowothers==1,
         suser='';
      else
         suser=[' AND tQueue.user="',c.owner,'"'];
      end
      %sql=['SELECT tQueue.*,tQ2.complete AS prevcomp',...
      %     ' FROM tQueue LEFT JOIN tQueue tQ2 ON tQueue.waitid=tQ2.id',...
      %     ' WHERE tQueue.complete<0 AND tQ2.complete=1'];
      sql=['SELECT -sum(tQueue.complete) as qactive,',...
           ' count(tQueue.complete) as qcount,tQueue.user',...
           ' FROM tQueue LEFT JOIN tQueue tQ2 ON tQueue.waitid=tQ2.id',...
           ' WHERE tQueue.complete in (-1,0)',...
           ' AND tQueue.allowqueuemaster=1',...
           ' AND (isnull(tQ2.complete) OR tQ2.complete=1)',...
           suser,...
           ' GROUP BY tQueue.user',...
           ' HAVING qcount-qactive>0',...
           ' ORDER BY qactive,qcount DESC'];
      userdata=mysql(sql);
      if length(userdata)==0,
         compok=0;
      end
      nextid=-1;
      
      if compok,
         % proceed with sending out a new job
         hostname=c.name; % [c.name,'.',c.ext];
         
         for ii=1:length(userdata),
            fprintf('user %s count %d/%d\n',userdata(ii).user, ...
                    ifstr2num(userdata(ii).qactive),userdata(ii).qcount);
         end
         
         uid=0;
         while nextid<0 & uid<length(userdata),
            uid=uid+1;
            gooduser=1;
            % ivar throttle - disabled 
            if 0 && strcmpi(userdata(uid).user,'ivar'),
                sql=['SELECT count(tQueue.id) as jobcount',...
                     ' FROM tQueue',...
                     ' WHERE complete=-1',...
                     ' AND computerid=',num2str(c.id),...
                     ' AND tQueue.user="',userdata(uid).user,'"'];
                udata=mysql(sql);
                if ~isempty(udata) && udata.jobcount>=2,
                    disp('skipping Ivar job temporarily');
                    gooduser=0;
                end
            end
            if gooduser,
                fprintf('checking for job for %s\n',userdata(uid).user);
                sql=['SELECT min(tQueue.id) as minid,',...
                     ' tQueue.user,tQueue.priority',...
                     ' FROM tQueue LEFT JOIN tQueue tQ2',...
                     ' ON tQueue.waitid=tQ2.id',...
                     ' WHERE tQueue.complete=0',...
                     ' AND tQueue.user="',userdata(uid).user,'"',...
                     ' AND tQueue.allowqueuemaster=1',...
                     ' AND (isnull(tQ2.complete) OR tQ2.complete=1)',...
                     ' GROUP BY tQueue.user,tQueue.priority',...
                     ' ORDER BY tQueue.priority DESC'];
                nextdata=mysql(sql);
                
                if length(nextdata)>0,
                    nextid=nextdata(1).minid;
                end
            end
         end
         
         if nextid>0,
            
            % ok, we're ready to go. make sure the machine is still alive
            l1=dbupdateload(hostname,'svd');
            if l1<0,
               dblog('host %s: pre-job check failed, marking dead',...
                     hostname);
               %dbset('tComputer',c.id,'allowqueuemaster',0);
               dbset('tComputer',c.id,'dead',1);
               dbevent(7,0,c.id);
            else
               
               % grab queue entry
               queuedata=dbgetqueue(nextid);
               r=dbsetqueue(nextid,0,-1,hostname);
               
               % make sure it was grabbed before someone else got it
               if r,
                  dbset('tQueue',0,'pid',nextid);
                  dbset('tComputer',compdata(compidx).id,'lastoverload',0);
                  
                  dblog('qid=%d: starting on %s user=%s progname=%s',...
                        queuedata.id,hostname,...
                        queuedata.user,queuedata.progname);
                  
                  if exist([LOGPATH,num2str(queuedata.id),'.out'],'file'),
                     delete([LOGPATH,num2str(queuedata.id),'.out']);
                  end
                  
                  shhostname=strsep(hostname,'.');
                  shhostname=shhostname{1};
                  % start the job on hostname via ssh
                  cmd=sprintf([SSHCMD,' -f ',queuedata.user,'@',shhostname,...
                              ' "',BINPATH,'runqueue ',...
                              num2str(queuedata.id),...
                              ' ', queuedata.progname,'"'''],...
                              queuedata.user);
                  
                  disp(cmd);
                  [s,w]=unix(cmd);
                  
                  % check to make sure it started correctly.
                  if s>0,
                     dblog('qid=%d: error starting. permissions problem?',...
                           queuedata.id);
                     r=dbsetqueue(nextid,0,2,hostname);
                  else
                     pause(0.01);
                     queuedata=dbgetqueue(nextid);
                     if queuedata.pid==0,
                        pause(0.1);
                        queuedata=dbgetqueue(nextid);
                     end
                     dbevent(2,queuedata.id,hostname,...
                             queuedata.pid,queuedata.user);
                  end
               end
            end
         end
      end
   end
   
   unix(''); % flush buffer
   
   %
   % fairness processing
   %
   
   if nextid<0,
      
      % cluster isn't full. fine.
      fprintf('not at queue capacity\n');
      pause(0.5);
      
   elseif mod(minute(now)-lastfairkillmin,60)>=10
      %
      % cluster is at capacity. figure out if anyone is hogging the
      % queue and kill one of their jobs if they are
      %
      sql=['SELECT user,sum(complete=-1) as runningjobs,',...
           ' sum(complete=0) as pendingdjobs',...
           ' FROM tQueue',...
           ' WHERE complete in (-1,0)',...
           ' GROUP BY user ',...
           ' ORDER BY runningjobs DESC'];
      jobdata=mysql(sql);
      
      actjobs=zeros(size(jobdata));
      pendjobs=actjobs;
      for jj=1:length(jobdata),
          actjobs(jj)=str2num(jobdata(jj).runningjobs);
          pendjobs(jj)=str2num(jobdata(jj).pendingdjobs);
      end
      totaljobs=sum([actjobs;0]);
      idealjobs=ceil(totaljobs./length(jobdata));

      cheatedusers=sum(pendjobs & actjobs<idealjobs);
      maxuser=min(find(actjobs==max(actjobs)));
      if actjobs(maxuser)>idealjobs+1 & cheatedusers>0,
         sql=['SELECT * FROM tQueue',...
              ' WHERE user="',jobdata(maxuser).user,'"',...
              ' AND complete=-1',...
              ' ORDER BY startdate DESC'];
         tq=mysql(sql);
         if length(tq)>1,
            dbset('tQueue',tq(1).id,'killnow',1);
            dblog('qid %d: user %s killing on %s to balance user load',...
                  tq(1).id,tq(1).user,tq(1).machinename);
            dbevent(5,tq(1).id,tq(1).machinename,tq(1).pid,tq(1).user);
         end
         lastfairkillmin=minute(now);
      end
      
      if totaljobs>0,
         sql=['SELECT user,count(id) as newjobs',...
              ' FROM tQueue WHERE complete=0',...
              ' AND user<>"',jobdata(1).user,'"',...
              ' GROUP BY user ',...
              ' ORDER BY newjobs DESC'];
         newjobdata=mysql(sql);
         
         newjobs=sum([cat(1,newjobdata.newjobs);0]);
         
         if totaljobs-jobdata(1).runningjobs<length(newjobs),
            
            % should kill a job so that everyone has at least one
            % running
            %disp('should kill a job');
         end
      end
   end
   
   % tell the database dbqueuemaster is still alive
   [s,dhost]=unix('hostname'); % flush
   [s,dhost]=unix('hostname');
   dhost=dhost(1:(end-1));
   sql=['UPDATE tGlobalData SET ',...
        'daemonclick=now(),daemonhost="',dhost,'"'];
   mysql(sql);
   
   % let things percolate
   pause(0.2);
   
   % periodically log current status
   if hour(now)~=lastreporthour | lastreporthour==-1,
      sql=['SELECT complete,count(id) as activecount',...
           ' FROM tQueue GROUP BY complete ORDER BY complete'];
      jobdata=mysql(sql);
      
      comp=cat(2,jobdata.complete);
      actcount=zeros(1,4);
      for cc=1:4,
         ccidx=find(comp==cc-2);
         if ~isempty(ccidx),
            actcount(cc)=jobdata(ccidx).activecount;
         end
      end
      
      dbevent(11,actcount(1));
      dbevent(12,actcount(2));
      dbevent(13,actcount(3));
      dbevent(14,actcount(4));
      
      sql=['SELECT count(load1) as nodecount,sum(dead) as deadcount,',...
           'sum(lastoverload * (1-dead)) as oloadcount,',...
           'sum(maxproc * (1-lastoverload) * (1-dead)) as maxproc,',...
           'avg(load1) as meanload FROM tComputer',...
           ' WHERE allowqueuemaster in (1,2);'];
      compdata=mysql(sql);
      
      dbevent(15,compdata.nodecount);
      dbevent(16,compdata.deadcount);
      dbevent(17,compdata.oloadcount);
      dbevent(18,round(compdata.meanload.*1000));
      
      fprintf(['act=%d,pend=%d,done=%d,dead=%d,',...
               'nodes=%d,dead=%d,olod=%d,load=%.2f\n'],...
              actcount,compdata.nodecount,compdata.deadcount,...
              compdata.oloadcount,compdata.meanload);
      lastreporthour=hour(now);
      
      calcheck;
   end
end


%catch
   
   % if queue crashes, email someone
   %cmd='fprintf(''Queue crashed!!!\n\n'');'
   %sub='Queue daemon crashed. Need to restart as user queued on nutmeg.';
   %emailres('prenger@socrates.berkeley.edu',cmd,sub);
   

%end
