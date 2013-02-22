% function dbcleanqueue(forceall)
%
% reset inactive, incomplete entries in tQueue.  also, always reset
% any -1's started by the matlab process that actually calls
% dbcleanqueue
%
% logic of cleaning:
% 1.  set all queue entries with complete=-1, pid=getpid, 
%                           and machinename=localhost to complete=0
% 2. for all queue entries with complete==-1,
%      if last update to tQueue was more than MAXLAT seconds ago,
%         if on the jlg subnet,
%            run dbupdate load to remotely check if process is still alive
%         otherwise
%            assume it's a dead job and set complete to 0
%
% optional parameter:
% forceall - reset all incomplete entries in tQueue.  ie, set all
%            -1's to 0.  BE CAREFUL!
%
% created SVD 6/2001
%
function dbcleanqueue(forceall)

dbopen;

MAXLAT=3*3600;  % number of seconds of inactivity to allow before
                % resetting the queue entry

% forceall effectively sets MAXLAT to 0
if not(exist('forceall')),
   forceall=0;
end

if forceall,
   disp('Force all entries to reset selected!');
   disp('Press a key to contiue or ctrl-C to cancel.');
   pause
end

% clear jobs for current matlab process
pid=getpid;

% figure out machinename: shortname + subnet (ie, jlg,bic,millennium)
machinename=getenv('MYHOST');
x=min([length(machinename)+1 findstr(machinename,'.')]);
shorthost=machinename(1:(x-1));
if x<length(machinename),
   ext=machinename((x+1):end);
   y=min([length(ext)+1 findstr(ext,'.')]);
   ext=ext(1:(y-1));
else
   ext='jlg';
end
machinename=[shorthost,'.',ext];

% force queue entries matching pid and machinename to be 0
sql=['UPDATE tQueue SET complete=0,progress=0,',...
           ' pid=0,machinename="",startdate=NULL',...
     ' WHERE pid=',num2str(pid),...
     ' AND machinename="',machinename,'"',...
     ' AND complete < 0'];
mysql(sql);

sql=['SELECT *,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
     ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
     ' FROM tQueue ORDER BY id'];
queuedata=mysql(sql);

for ii=1:length(queuedata),
   
   % check if queue progress has been updated recently
   if (queuedata(ii).secago>MAXLAT | forceall) & queuedata(ii).complete==-1,
      if findstr(queuedata(ii).machinename,'.jlg'),
         % if it's on the local subnet, we can check to see if the
         % process is still alive via rsh
         dbupdateload(queuedata(ii).machinename);
      else
         % otherwise, we just have to assume it's dead
         fprintf('Queue id %d inactive %d sec.\n', ...
                 queuedata(ii).id,queuedata(ii).secago);
         dbsetqueue(queuedata(ii).id,0,0);
      end
   else
      %fprintf('queue id %d ok (%d sec).\n', ...
      %        queuedata(ii).id,queuedata(ii).secago);
   end
end

   


