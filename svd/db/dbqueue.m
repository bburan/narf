% function qid=dbqueue([complete=<all>],notemask[='%'],[user=getenv('USER')])
%
% list entries in tQueue analysis queue
%
% parameters:
% complete - filter by completion status
%            -1 = in progress
%             0 = not yet run
%             1 = complete
%             2 = dead
%      'string' = active jobs with note or user containing string
% user - list only entries matching user.  generally, an individual
%        only wants to list their own jobs.  unix login id is
%        stored in the global DBUSER upon excution of dbopen().
%        user='%' will return all entries
% notemask - list only entries with a note field "like" notemask.
%            notemask='%' will return all entries
%
% created SVD 6/2001
%
function qid=dbqueue(complete,notemask,user)

global DBUSER

dbopen;

if ~exist('notemask','var'),
   notemask='%';
end
if ~exist('user','var') | isempty(user),
   user=DBUSER;
end

if exist('complete','var') & ~isnumeric(complete),
   if strcmp(complete,'all'),
      user='';
      notemask='';
   else
      user=complete;
      notemask=complete;
   end
   sql=['SELECT *,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
        ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
        ' FROM tQueue',...
        ' WHERE complete=-1',...
        ' AND (user like "%',user,'%"'...
        ' OR note like "%',notemask,'%")'...
        ' ORDER BY id'];
elseif exist('complete','var'),
   sql=['SELECT *,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
        ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
        ' FROM tQueue',...
        ' WHERE complete=',num2str(complete),...
        ' AND user like "%',user,'%"'...
        ' AND note like "%',notemask,'%"'...
        ' ORDER BY id'];
else
   sql=['SELECT *,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
        ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
        ' FROM tQueue',...
        ' WHERE user like "%',user,'%"'...
        ' AND note like "%',notemask,'%"'...
        ' ORDER BY id'];
end
queuedata=mysql(sql);

if length(queuedata)>0,
   progname=queuedata(1).progname;
else
   progname='';
end
if strcmp(progname,'free2mov') | strcmp(progname,'cellxc') | ...
      strcmp(progname,'kerncomp'),
   fprintf('   ID. CELL  (BA) RUN  CP    PID ROUTINE     USER     PROG LAST MACHINE\n');
else
   fprintf('   ID. CP   PID NOTES       USER    COUNT LAST MACHINE\n');   
end

for ii=1:length(queuedata),
   
   if queuedata(ii).complete==1,
      machinename='--';
      pid=0;
      progress=queuedata(ii).progress;
      secago=0;
   elseif queuedata(ii).complete==-1,
      machinename=queuedata(ii).machinename;
      pid=queuedata(ii).pid;
      progress=queuedata(ii).progress;
      secago=queuedata(ii).secago;
   elseif queuedata(ii).complete==2,
      machinename='-- *DEAD*';
      pid=0;
      progress=0;
      secago=0;
   else
      machinename='--';
      pid=0;
      progress=0;
      secago=0;
   end
   
   if queuedata(ii).killnow,
      sextra=' *KILL PENDING*';
   else
      sextra='';
   end
   
   progname=queuedata(ii).progname(1:min([11 length(queuedata(ii).progname)]));
   if ~isempty(queuedata(ii).note),
      note=queuedata(ii).note;
   else
      note=progname;
   end
   
   fprintf('%4d. %2d %6d %-11s %-8s%5d %4d %s%s\n',...
           queuedata(ii).id,queuedata(ii).complete,pid,...
           note,queuedata(ii).user,progress,secago,machinename,sextra);
end

if nargout>0,
   qid=cat(1,queuedata.id);
end



