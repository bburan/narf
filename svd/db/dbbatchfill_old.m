% function runids=dbbatchfill(batchid,matchbatchid,yestoall,sarea);
function runids=dbbatchfill(batchid,matchbatchid,yestoall,sarea);

MINLEN=0; % minimum total length of respfiles allowed into
          % batch. 0 lets everyone in.
          
if ~exist('yestoall','var'),
   yestoall=0;
end

dbopen;

sql=['SELECT * FROM sBatch WHERE id=',num2str(batchid)];
params=mysql(sql);
if length(params.parmstring)>0,
   eval(char(params.parmstring));
end

rejectcells='("R240A","R247B","modelss")';

if 1 || params.stimspeedid==0,
   disp('disabling speed test');
   speedtest='';
elseif params.stimspeedid>=60,
   speedtest=[' AND stimspeedid>=',num2str(params.stimspeedid)];
else
   speedtest=[' AND stimspeedid=',num2str(params.stimspeedid)];
end

if ~isfield(params,'stimsnr') || isempty(params.stimsnr),
   snrtest='';
elseif params.stimsnr>=100,
   snrtest=[' AND stimsnr>=',num2str(params.stimsnr)]
else
   snrtest=[' AND stimsnr=',num2str(params.stimsnr)];
end
if isfield(params,'area'),
   fprintf('restricted to cells matching area="%s"\n',params.area);
   sarea=[' AND sCellFile.area like "',params.area,'"'];
else
   sarea='';
end

switch batchid,
 case 139,
  disp('batch 139: special, only including long CCHs');
  speedtest=[' AND stimspeedid>=1'];
 
 case 143,
  disp('batch 143: special, only including first runclassid');
  params.runclassid=params.runclassid(1);
end

rcstr='runclassid in (';
for ii=1:length(params.runclassid),
   rcstr=[rcstr,num2str(params.runclassid(ii)),','];
end
rcstr(end)=')';

%keyboard
if exist('matchbatchid','var') && ~isempty(matchbatchid),
   fprintf('restricted to cells in batch %d\n',matchbatchid);
   sql=['SELECT DISTINCT sCellFile.singleid,sCellFile.masterid,',...
        ' sCellFile.cellid,sCellFile.path,sum(resplen) as totlen',...
        ' FROM sCellFile INNER JOIN gSingleRaw',...
        ' ON sCellFile.singlerawid=gSingleRaw.id',...
        ' INNER JOIN sRunData ON (sRunData.batch=',num2str(matchbatchid),...
        ' AND sRunData.cellid=sCellFile.cellid)',...
        ' WHERE ',rcstr,speedtest,snrtest,sarea,...
        ' AND respfmtcode=',num2str(params.respfmtcode),...
        ' AND stimfmtcode=',num2str(params.stimfmtcode),...
        ' AND not(sCellFile.cellid in ',rejectcells,')',...
        ' AND not(gSingleRaw.crap)',...
        ' GROUP BY sCellFile.cellid',...
        ' ORDER BY sCellFile.cellid,path'];
else
   area='';
   sql=['SELECT DISTINCT sCellFile.singleid,sCellFile.masterid,',...
        ' sCellFile.cellid,sCellFile.path,sum(resplen) as totlen',...
        ' FROM sCellFile INNER JOIN gSingleRaw',...
        ' ON sCellFile.singlerawid=gSingleRaw.id',...
        ' WHERE ',rcstr,speedtest,snrtest,sarea,...
        ' AND respfmtcode=',num2str(params.respfmtcode),...
        ' AND stimfmtcode=',num2str(params.stimfmtcode),...
        ' AND not(sCellFile.cellid in ',rejectcells,')',...
        ' AND not(gSingleRaw.crap)',...
        ' GROUP BY sCellFile.cellid',...
        ' ORDER BY sCellFile.cellid,path'];
end
filedata=mysql(sql);

% find current entries in sRunData that don't have any sCellFiles
sql=['SELECT DISTINCT sRunData.* FROM sRunData LEFT JOIN sCellFile',...
     ' ON sRunData.singleid=sCellFile.singleid',...
     ' WHERE ',rcstr,speedtest,...
     ' AND respfmtcode=',num2str(params.respfmtcode),...
     ' AND stimfmtcode=',num2str(params.stimfmtcode),...
     ' AND isnull(sCellFile.id)',...
     ' AND batch=',num2str(batchid)];
badrundata=mysql(sql);

rcsetstrings;
if ~isempty(params.expfrac) & params.expfrac>0 & params.expfrac<1,
   sspeedfrac=sprintf('.%0.3d',round(params.expfrac*100));
elseif ~isempty(params.stimspeedid) & params.stimspeedid>0,
   sspeedfrac=['.',num2str(params.stimspeedid)];
else
   sspeedfrac='';
end

% switch from jlg to nsl labs
if batchid<108,
   respath=['/auto/k5/david/data/batch',num2str(batchid),'/'];
else
   respath=['/auto/data/nsl/users/svd/data/batch',...
            num2str(batchid),'/'];
end

samerespath=1;
sql=['SELECT * FROM sRunData WHERE batch=',num2str(batchid)];
rundata=mysql(sql);
if length(rundata)>=2,
   if ~strcmp(rundata(1).respath,rundata(2).respath),
      samerespath=0;
      respath='*scellfile-path*';
   end
end

for ii=1:length(rundata),
   sql=['SELECT sCellFile.*,gSingleRaw.crap',...
        ' FROM sCellFile INNER JOIN gSingleRaw',...
        ' ON sCellFile.singlerawid=gSingleRaw.id',...
        ' WHERE ',rcstr,speedtest,...
        ' AND respfmtcode=',num2str(params.respfmtcode),...
        ' AND stimfmtcode=',num2str(params.stimfmtcode),...
        ' AND sCellFile.cellid="',rundata(ii).cellid,'"'];
   scdata=mysql(sql);
   if length(scdata)==0,
      fprintf('cell %s missing... deleting sRunData entry\n',...
              rundata(ii).cellid);
      sql=['DELETE FROM sRunData WHERE id=', ...
           num2str(rundata(ii).id)];
      mysql(sql);
      sql=['DELETE FROM sResults WHERE runid=', ...
           num2str(rundata(ii).id)];
      mysql(sql);
   elseif sum(1-cat(1,scdata.crap))==0,
      fprintf('file for cell %s is crap... deleting sRunData entry\n',...
              rundata(ii).cellid);
      sql=['DELETE FROM sRunData WHERE id=', ...
           num2str(rundata(ii).id)];
      mysql(sql);
      sql=['DELETE FROM sResults WHERE runid=', ...
           num2str(rundata(ii).id)];
      mysql(sql);
   end
end

sql=['SELECT * FROM gRunClass WHERE id=',num2str(params.runclassid(1))];
runclassdata=mysql(sql);
runclass=runclassdata.name;

resbase=['.',runclass,'.',params.kernfmt, ...
         '.',num2str(batchid),'.res.mat'];
kernbase=['.',runclass,'.',params.kernfmt, ...
          '.',num2str(batchid),'.kern.mat'];

fprintf('respath: %s\nresbase: %s\nkernbase: %s\n',...
        respath,resbase,kernbase);
%disp('paused before adding new runs');
%keyboard

runids=[];
runcount=0;
for ii=1:length(filedata),
   if filedata(ii).totlen>MINLEN,
      sql=['SELECT * FROM sRunData',...
           ' WHERE singleid=',num2str(filedata(ii).singleid),...
           ' AND batch=',num2str(batchid)];
      rundata=mysql(sql);
      
      if length(rundata)>0,
         %fprintf('Entry exists for cellid %s. Skipping.\n',...
         %        filedata(ii).cellid);
         
      else
         if yestoall,
            fprintf('Auto-adding entry for cellid %s\n',filedata(ii).cellid);
            yn='y';
         else
            yn=input(sprintf('Add entry for cellid %s (y/[n]) ',...
                             filedata(ii).cellid),'s');
         end
         if ~isempty(yn) & strcmp(yn(1),'y'),
            if ~samerespath,
               respath=filedata(ii).path;
            end
            
            sql=['INSERT INTO sRunData ',...
                 ' (singleid,masterid,cellid,batch,respath,resfile,kernfile)',...
                 ' VALUES (',num2str(filedata(ii).singleid),',',...
                 num2str(filedata(ii).masterid),',',...
                 '"',filedata(ii).cellid,'",',...
                 num2str(batchid),',',...
                 '"',respath,'",',...
                 '"',[filedata(ii).cellid,resbase],'",',...
                 '"',[filedata(ii).cellid,kernbase],'")'];
            [aa,bb,newrunid]=mysql(sql);
            runids=[runids;newrunid];
            runcount=runcount+1;
         end
      end
   end
end

fprintf('batchid=%d added %d runs, resbase=%s\n',...
        batchid,runcount,resbase);


