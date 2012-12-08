% function r=dbdeletequeue(queueidx,[user=getenv('USER')]);
%
% remove entry with id=queueidx from tQueue
% 
% user parameter must match user in db entry
% r=1 if deleted successfully, r=0 if not.
%
% created SVD 1/23/2003
%
function r=dbdeletequeue(queueidx,user);

global DBUSER
if ~exist('user','var'),
   user=DBUSER;
end

dbopen;

r=0;
queuedata=dbgetqueue(queueidx);

if length(queuedata)>0 & (strcmp(user,'%') | strcmp(user,queuedata.user)),
   
   % user matches. proceed.
   if queuedata.complete==-1,
      fprintf('Killing and deleting queue id %d.\n',queueidx);
      dbkillqueue(queuedata.id,3);
      pause(0.5);
   else
      fprintf('Deleting queue id %d.\n',queueidx);
      sql=['DELETE FROM tQueue WHERE id=',num2str(queueidx)];
      mysql(sql);
   end
   r=1;
   
else
   disp('no match');
end






