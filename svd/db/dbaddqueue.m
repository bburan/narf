% function queueidx=dbaddqueue(rundataid,progname,parmstring,[allowqueuemaster],[note],[user]);
%
% add an entry to tQueue
%
% Created SVD 6/01
%
function queueidx=dbaddqueue(rundataid,progname,parmstring,allowqueuemaster,note,user);

dbopen;

global DBUSER
if ~exist('allowqueuemaster','var'),
   allowqueuemaster=0;
end
if ~exist('note','var'),
   note='';
end
if ~exist('user','var'),
   user=DBUSER;
end

y=year(now);
mon=month(now);
d=day(now);
h=hour(now);
min=minute(now);
s=round(second(now));
queuedate=sprintf('%d-%.2d-%.2d %.2d:%.2d:%.2d',y,mon,d,h,min,s);

[aff,queueidx]=sqlinsert('tQueue',...
                         'rundataid',rundataid,...
                         'progname',progname,...
                         'parmstring',parmstring,...
                         'priority',1,...
                         'queuedate',queuedate,...
                         'allowqueuemaster',allowqueuemaster,...
                         'note',note,...
                         'user',user);
   
fprintf('Created queue entry %d.\n',queueidx);
dbevent(6,queueidx);
