% function id=dbevent(code,queueid,hostname/computerid,pid,userid,note);
%
% event codes:
% 1 - queue daemon start/status
% 2 - job start
% 3 - job finish
% 4 - job die
% 5 - job kill
% 6 - job add
% 7 - machine die
%
% hourly status update "events" -- values stored in queueid
% 11 - active job count
% 12 - pending job count
% 13 - completed jobs
% 14 - dead jobs
% 15 - node count
% 16 - dead nodes
% 17 - overloaded nodes
% 18 - mean load (*1000) 
%
%
function id=dbevent(code,queueid,computerid,pid,userid,note);

if ~exist('code','var'),
   disp('must specify event code');
   return
end
if ~exist('queueid','var'),
   global BATQUEUEID
   if ~isempty(BATQUEUEID),
      queueid=BATQUEUEID;
   else
      queueid=0;
   end
end
if ~exist('computerid','var'),
   computerid=dbgetcompid;
end
if ~isnumeric(computerid),
   computerid=dbgetcompid(computerid);
end
if ~exist('pid','var'),
   pid=getpid;
end
if ~exist('userid','var'),
   userid=getenv('USER');
end

if ~exist('note','var'),
   [aff,id]=sqlinsert('tEvent','code',code,'queueid',queueid,...
                      'computerid',computerid,'pid',pid,'user',userid);
else
   [aff,id]=sqlinsert('tEvent','code',code,'queueid',queueid,...
                      'computerid',computerid,'pid',pid,'user',userid,...
                      'note',note);
end






