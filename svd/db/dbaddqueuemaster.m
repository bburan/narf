% function queueidx=dbaddqueuemaster(parmstring,[note],priority[=1]);
%
% add a job to the queue
%
% parameters: 
% parmstring - string containing matlab command(s). executed
%              exactly as if pasted at the matlab prompt
%              immediately after start-up in your home directory.
%
%              dependent sequences: if parmstring is a cell array
%              of strings (eg, parmstring={'prog1;','prog2;'}) a
%              job will be added for each array entry. a given job
%              will not start until the previous job has completed
%              sucessfully.
%              
% note - note to appear for this job on the queue monitor
% priority - (0-5) higher values cause jobs to be executed sooner than
%            jobs with lower values (among your OWN jobs)
%
% created SVD 6/7/03
% modified SVD 7/28/04 - added dependent job sequence options.
%
function queueidx=dbaddqueuemaster(parmstring,note,priority,rundataid,user);

global DBISOPEN DBUSER

dbopen;

if ~exist('user','var'),
   user=DBUSER;
   if isempty(user),
      user=getenv('USER');
   end
   if isempty(user),
      user=getenv('user');
   end
end
if ~exist('note','var'),
   note='';
end
if ~exist('priority','var'),
   priority=1;
end
if ~exist('rundataid','var'),
   rundataid=0;
end

progname='matlabbg queuerun';
allowqueuemaster=1;

if iscell(parmstring),
   seqcount=length(parmstring);
else
   seqcount=1;
   parmstring={parmstring};
end

y=year(now);
mon=month(now);
d=day(now);
h=hour(now);
min=minute(now);
s=round(second(now));
queuedate=sprintf('%d-%.2d-%.2d %.2d:%.2d:%.2d',y,mon,d,h,min,s);

waitid=0;
for seqidx=1:seqcount,
   if iscell(note),
      tnote=note{seqidx};
   elseif isempty(note) | seqcount==1,
      tnote=note;
   else
      tnote=sprintf('%s(%d)',note,seqidx);
   end
   
   [aff,queueidx]=sqlinsert('tQueue',...
                            'rundataid',rundataid,...
                            'progname',progname,...
                            'priority',priority,...
                            'parmstring',parmstring{seqidx},...
                            'queuedate',queuedate,...
                            'allowqueuemaster',allowqueuemaster,...
                            'user',user,...
                            'note',tnote,...
                            'waitid',waitid);
   waitid=queueidx;
   fprintf('Created queue entry %d.\n',queueidx);
   dbevent(6,queueidx);
end

