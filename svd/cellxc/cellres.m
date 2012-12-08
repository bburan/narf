% function r=cellres(cellid,batch,savefigs)
%
% load results from resfile in sRunData and display fits, pred
% results
%
% r=0 if no entries found in db, =1 otherwise
%
function r=cellres(runidx,batch,savefigs)

if ~exist('savefigs','var'),
    savefigs=0;
end

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
        ' AND batch=',num2str(batch)];
   rundata=mysql(sql);
end

if length(rundata)==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

xcresult([rundata(1).respath,rundata(1).resfile],[],savefigs);


% special case for TORCs:
if batch==108,
   rctorc(cellid,batch);
end