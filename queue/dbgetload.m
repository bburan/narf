% function [r,compnames]=dbgetload(computerid,silent,updatethresh,updateloadbg)
%
% r is five columns: (id, load1, numproc, location, secsinceupdate)
%
% optional parameters:
%
function [r,compnames]=dbgetload(computerid,silent,updatethresh,updateloadbg)

if ~exist('updatethresh','var')
   updatethresh=30;    % update load every updatethresh seconds.
end
if ~exist('updateloadbg','var'),
   updateloadbg=0;
end

dbopen;
if ~exist('silent','var'),
   silent=0;
end
if exist('computerid','var') & ~isempty(computerid),
   sql=['SELECT *,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
        ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
        ' FROM tComputer where id=',num2str(computerid),...
        ' order by load1,name'];
else
   sql=['SELECT * ,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
       ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
       ' FROM tComputer order by load1,name'];
end
compdata=mysql(sql);

r=zeros(length(compdata),5);
compnames={};
for ii=1:length(compdata),
   
   compnames{ii}=compdata(ii).name;
   r(ii,1)=compdata(ii).id;
   r(ii,2)=compdata(ii).load1;
   r(ii,3)=compdata(ii).numproc;
   r(ii,4)=compdata(ii).location;
   r(ii,5)=compdata(ii).secago;
   %r(ii,6)=compdata(ii).allowqueuemaster;
   
   if compdata(ii).numproc==0,
      tupdatethresh=updatethresh*4;
   else
      tupdatethresh=updatethresh;
   end
   
   % changed to allowqueuemaster >= 0 to include "off" computers in
   % load update checks
   if compdata(ii).allowqueuemaster>=0 & ...
          compdata(ii).nocheck==0 & ...
          compdata(ii).dead==0 & ...
          compdata(ii).secago>tupdatethresh,
      if updateloadbg
         fprintf('.');
      else
         fprintf([compnames{ii},'.']);
      end
      dbupdateload([compnames{ii},'.',compdata(ii).ext],...
                   getenv('USER'),updateloadbg);
   end
   
   if compdata(ii).secago >= 10000,
      % don't bother displaying
      %fprintf('%5d ',99999);
      r(ii,3)=0;
   else
      if ~silent,
         fprintf('(%3d) %-17s %2d/%d %.2f\n',...
              r(ii,1),[compdata(ii).name,'.',compdata(ii).ext],...
              compdata(ii).numproc,compdata(ii).maxproc,r(ii,2));
      end
   end
end






