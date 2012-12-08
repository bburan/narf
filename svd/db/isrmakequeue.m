% function isrmakequeue(cellid,batch,progname,pausetime)
%
% quick job adder for isr queue%
% 
% created svd 2007-05-30 -- ripped off cellxcmakequeue.m
%
function isrmakequeue(cellid,batch,progname,pausetime);

dbopen;

if not(exist('cellid','var')), % do individual run? if > 0
   error('must specify cellid');
end
if ~isnumeric(cellid) && not(exist('batch','var')),
   error('must specify batch');
end

if ~exist('progname','var'),
   if ~isnumeric(cellid),
      batchdata=dbget('sBatch',batch);
      progname=batchdata.matcmd;
   else
      progname='';
   end
end
if ~exist('pausetime','var')
   pausetime=0;
end

allowqueuemaster=0;
parmstring='';

if ~isempty(cellid) && ~isnumeric(cellid),
   sql=['SELECT * FROM sRunData where cellid="',cellid,'"',...
        ' AND batch=',num2str(batch)];
elseif ~isempty(cellid) && isnumeric(cellid),
   rstr=mat2str(cellid(:)');
   rstr(1)='(';
   rstr(end)=')';
   rstr=strrep(rstr,' ',',');
   
   sql=['SELECT DISTINCT sRunData.*,matcmd FROM sRunData,sBatch',...
        ' WHERE sRunData.batch=sBatch.id AND sRunData.id in ',rstr];
   
else
   sql=['SELECT DISTINCT * FROM sRunData',...
        ' WHERE batch=',num2str(batch),...
        ' ORDER BY cellid,id'];
   %     ' AND not(cellid like "%model")',...
   % allow model cells in
end
rundata=mysql(sql);
runcount=length(rundata);

for rr=1:runcount,
   
   runidx=rundata(rr).id;
   % don't need this any more??
   %dbcleanqueue;
   
   note=sprintf('%s/%d',rundata(rr).cellid,rundata(rr).batch);
   
   if isempty(progname),
      tprogname=rundata(rr).matcmd;
   else
      tprogname=progname;
   end
   
   parmstring=sprintf('%s(''%s'',%d);',...
                      tprogname,rundata(rr).cellid,rundata(rr).batch);
   queueidx=dbaddqueueisr(parmstring,note);
   
   fprintf('Cell: %s  queueidx=%d\n',rundata(rr).cellid,queueidx);
   
   if pausetime>0,
      fprintf('pausing %d sec...\n',pausetime);
      pause(pausetime);
   end
   
end

fprintf('Added %d entries to queue.\n',runcount);
