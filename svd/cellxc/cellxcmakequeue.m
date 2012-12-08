% function cellxcmakequeue(runidx,rungroup,batch,progname)
%
% add entries to queue for movmatchsvd routine.  decorr paramters are
% hard-coded into makequeue.m, so you have to edit the m-file
% directly at this point.
%
% 2 cellxcmaster
% 3 kerncompmaster
% 4 xcfittune
% 6 xcstsep
% 7 xcsepfull
% 8 *RESERVED* -- make rungroup a string
%    eg, 'xcsepfull', 'xcstsep', 'kerncompmaster'
%
function cellxcmakequeue(runidx,rungroup,batch,progname);

dbopen;

if not(exist('runidx','var')), % do individual run? if > 0
   runidx=[];
end
if not(exist('rungroup','var')),
   rungroup=1;
end
if ~isnumeric(rungroup),
   matcmd=rungroup;
   rungroup=8;
end

if not(exist('batch','var')),
   batch=1;         % do whole batch if length(runidx)>0
end

if rungroup==1,
   progname='cellxc';
   allowqueuemaster=0;
   parmstring='';
elseif ismember(rungroup,[2 3 4 6 7 8]),
   progname='matlabbg queuerun';
   allowqueuemaster=1;
   parmstring='';
elseif rungroup==5,
   progname='matlabbg tnl';
   allowqueuemaster=1;
   parmstring='';
else
   if not(exist('progname','var')),
      progname='cellxc';
   end
   allowqueuemaster=0;
   parmstring='';
end

if length(runidx) > 0,
   sql=['SELECT * FROM sRunData where id in ('];
   for ii=1:length(runidx),
      if ii==1,
         sql=[sql,num2str(runidx(ii))];
      else
         sql=[sql,',',num2str(runidx(ii))];
      end
   end
   sql=[sql,');'];
else
   sql=['SELECT DISTINCT sRunData.*',...
        ' FROM sRunData',...
        ' WHERE sRunData.batch=',num2str(batch),...
        ' ORDER BY cellid,id'];
   %     ' AND not(cellid like "%model")',...
   % allow model cells in
end
rundata=mysql(sql);
runcount=length(rundata);

for rr=1:runcount,
   
   runidx=rundata(rr).id;
   dbcleanqueue;
   
   note=sprintf('%s/%d',rundata(rr).cellid,rundata(rr).batch);
   
   if ismember(rungroup,[1 5]),
      % old way of doing this
      queueidx=dbaddqueue(runidx,progname,parmstring,allowqueuemaster,note);
   
   elseif rungroup==2,
      parmstring=sprintf('cellxcmaster(''%s'',%d);',...
                         rundata(rr).cellid,rundata(rr).batch);
      queueidx=dbaddqueueisr(parmstring,note);
      
   elseif rungroup==3,
      parmstring=sprintf('kerncompmaster(''%s'',%d);',...
                         rundata(rr).cellid,rundata(rr).batch);
      queueidx=dbaddqueueisr(parmstring,note);
      
   elseif rungroup==4,
      parmstring=sprintf('xcfittune(''%s'',%d);',...
                         rundata(rr).cellid,rundata(rr).batch);
      queueidx=dbaddqueueisr(parmstring,note);
      
  elseif rungroup==6,
      parmstring=sprintf('xcstsep(''%s'',%d);',...
                         rundata(rr).cellid,rundata(rr).batch);
      queueidx=dbaddqueueisr(parmstring,note);
      
   elseif rungroup==7,
      parmstring=sprintf('xcsepfull(''%s'',%d);',...
                         rundata(rr).cellid,rundata(rr).batch);
      queueidx=dbaddqueueisr(parmstring,note);
      
   elseif rungroup==8,
      % new generic way of doing it. much easier.
      parmstring=sprintf('%s(''%s'',%d);',...
                         matcmd,rundata(rr).cellid,rundata(rr).batch);
      queueidx=dbaddqueueisr(parmstring,note);
      
   else
      disp('rungroup not defined!!!');
      queueidx=-1;
   end
   fprintf('Cell: %s  queueidx=%d\n',rundata(rr).cellid,queueidx);   
end

if runcount>1 & ismember(batch,[23 24 26 27 28 29]),
   mailto='svd@umd.edu';
   mailcommand=['dbresv1(',num2str(batch),')'];
   sql=['UPDATE tQueue set',...
        ' mailto="',mailto,'",',...
        ' mailcommand="',mailcommand,'"',...
        ' WHERE id=',num2str(queueidx)];
   mysql(sql);
end

fprintf('Added %d entries to queue.\n',runcount);
