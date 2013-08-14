% function r=dbtickqueue(extended_status);
%
% tick progress counter in tQueue entry, optionally save an
% extended_status double(). automatically set complete=-1 as well
%
% created SVD 8/2013, ripped off of dbsetqueue
%
function r=dbtickqueue(extended_status);

global BATQUEUEID
global USEDB

if ~USEDB,
   r=1;
   return;
end

persistent lasttic

if now-lasttic < 1e-4,
    r=1;
    %now-lasttic
    return
end

lasttic=now;

if isempty(BATQUEUEID) | BATQUEUEID<=0,
   return
end
if ~exist('extended_status') || isempty(extended_status),
    extended_status=0;
end

dbopen;

sql=['UPDATE tQueue ',...
     'SET progress=progress+1,',...
     ' extended_status=',num2str(extended_status),',',...
     ' complete=-1',...
     ' WHERE id=',num2str(BATQUEUEID)];

[res,r]=mysql(sql);

% always update load to keep dbqueuemaster happy
dbupdateload;
