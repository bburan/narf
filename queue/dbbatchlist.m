%function res=dbbatchlist(batchnum)
%
function res=dbbatchlist(batchnum)

dbopen;

if exist('batchnum','var') & isnumeric(batchnum),
   sql=['SELECT * FROM sBatch where id=',num2str(batchnum)];
   params=mysql(sql);
   if length(params.parmstring)>0,
      eval(char(params.parmstring));
   end
   
   if isempty(params.stimspeedid) || params.stimspeedid==0,
      speedtest='';
   elseif params.stimspeedid>=60,
      speedtest=[' AND stimspeedid>=',num2str(params.stimspeedid)];
   else
      speedtest=[' AND stimspeedid=',num2str(params.stimspeedid)];
   end
   if length(params.runclassid)>1,
      rcstr='sCellFile.runclassid in (';
      for ii=1:length(params.runclassid),
         rcstr=[rcstr,num2str(params.runclassid(ii)),','];
      end
      rcstr(end)=')';
   else
      rcstr=['sCellFile.runclassid=',num2str(params.runclassid)];
   end
   
   runclass=dbget('gRunClass',params.runclassid(1),'name');
   sql=['SELECT DISTINCT stimfilefmt FROM sCellFile where stimfmtcode=',...
        num2str(params.stimfmtcode)];
   data=mysql(sql);
   stimfilefmt=data.stimfilefmt;
   sql=['SELECT DISTINCT respfilefmt FROM sCellFile where respfmtcode=',...
        num2str(params.respfmtcode)];
   data=mysql(sql);
   respfilefmt=data(1).respfilefmt;
   
   sql=['SELECT sRunData.*,sum(sCellFile.resplen) as resplen',...
        ' FROM sRunData INNER JOIN sCellFile',...
        ' ON sRunData.singleid=sCellFile.singleid',...
        ' WHERE ',rcstr,...
        speedtest,...
        ' AND sCellFile.stimfmtcode=',num2str(params.stimfmtcode),...
        ' AND sCellFile.respfmtcode=',num2str(params.respfmtcode),...
        ' AND batch=',num2str(batchnum),...
        ' GROUP BY singleid',...
        ' ORDER BY sRunData.cellid']
   rundata=mysql(sql);
   
   fprintf('Batch %d: %s (%s)\n',batchnum,params.name,params.details);
   fprintf('Run class: %s',runclass);
   if params.stimspeedid>0,
      fprintf(' %d Hz',params.stimspeedid);
   end
   fprintf('\nStim fmt: %s Resp fmt: %s\n',stimfilefmt,respfilefmt);
   
   %keyboard
   
   fprintf('%-4s %-10s %-8s %-24s\n','ID','CELL','RESLEN','RESFILE');
   for runidx=1:length(rundata),
      
      fprintf('%-4d %-10s %8d %-24s\n',rundata(runidx).id,...
              rundata(runidx).cellid,ifstr2num(rundata(runidx).resplen),...
              rundata(runidx).resfile);
   end
   if nargout>0,
      res=rundata;
   end
else
   if ~exist('batchnum','var'),
      batchnum='';
   end
   
   sql=['SELECT sBatch.*, count(sRunData.id) as runcount ',...
        ' FROM sBatch INNER JOIN sRunData',...
        ' ON sBatch.id=sRunData.batch',...
        ' WHERE sBatch.name like "%',batchnum,'%"',...
        ' GROUP BY sBatch.id',...
        ' ORDER BY sBatch.id'];
   params=mysql(sql);
   
   fprintf('%-3s %-24s %3s %-15s (RN NW DN DD / TO)\n','ID','NAME','CT','CMD');
   for bidx=1:length(params),
      sql=['SELECT sum(complete=-1) as runcount,'...
           ' sum(complete=0) as newcount,',...
           ' sum(complete=1) as donecount,',...
           ' sum(complete=2) as deadcount,',...
           ' count(id) as totalcount',...
           ' FROM tQueue',...
           ' WHERE substring_index(note,"/",-1)="',...
           num2str(params(bidx).id),'"'];
      qdata=mysql(sql);
      if qdata.totalcount==0,
         fprintf('%-3d %-24s %3d %s\n',...
                 params(bidx).id,params(bidx).name,...
                 params(bidx).runcount,params(bidx).matcmd);
      else
         fprintf('%-3d %-24s %3d %-15s (%d %d %d %d / %d)\n',...
                 params(bidx).id,params(bidx).name,...
                 params(bidx).runcount,params(bidx).matcmd,...
                 str2num(qdata.runcount),str2num(qdata.newcount),...
                 str2num(qdata.donecount),str2num(qdata.deadcount),...
                 qdata.totalcount);
      end
   end
   if nargout>0,
      res=params;
   end
end

