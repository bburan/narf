% function runids=dbadd2batch(cellid,batchid);
function runid=dbadd2batch(cellid,batchid);

dbopen;

sql=['SELECT * FROM sRunData WHERE batch=',num2str(batchid),...
     ' AND cellid="',cellid,'"'];
rundata=mysql(sql);
if ~isempty(rundata),
   warning('entry already exists');
   return
end

sql=['SELECT * FROM sBatch WHERE id=',num2str(batchid)];
params=mysql(sql);
if length(params.parmstring)>0,
   eval(char(params.parmstring));
end

if ~isfield(params,'stimsnr') || isempty(params.stimsnr),
   snrtest='';
elseif params.stimsnr>=100,
   snrtest=[' AND stimsnr>=',num2str(params.stimsnr)]
else
   snrtest=[' AND stimsnr=',num2str(params.stimsnr)];
end
if isfield(params,'area'),
   sarea=params.area;
end

rcstr='runclassid in (';
for ii=1:length(params.runclassid),
   rcstr=[rcstr,num2str(params.runclassid(ii)),','];
end
rcstr(end)=')';

sarea='%';
fprintf('restricted to cells matching area="%s"\n',sarea);
sql=['SELECT DISTINCT sCellFile.singleid,sCellFile.masterid,',...
     ' sCellFile.cellid,sCellFile.path,',...
     ' sum(sCellFile.resplen) as totlen',...
     ' FROM sCellFile INNER JOIN gSingleRaw',...
     ' ON sCellFile.singlerawid=gSingleRaw.id',...
     ' WHERE ',rcstr,snrtest...
     ' AND respfmtcode=',num2str(params.respfmtcode),...
     ' AND stimfmtcode=',num2str(params.stimfmtcode),...
     ' AND not(gSingleRaw.crap)',...
     ' AND sCellFile.area like "',sarea,'"',...
     ' AND sCellFile.cellid="',cellid,'"',...
     ' GROUP BY sCellFile.cellid'];
filedata=mysql(sql);

if isempty(filedata),
   warning('no matches in sCellFile');
   return
end

% switch from jlg to nsl labs
if batchid<108,
   respath=['/auto/k5/david/data/batch',num2str(batchid),'/'];
else
   respath=['/auto/data/nsl/users/svd/data/batch',...
            num2str(batchid),'/'];
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

sql=['INSERT INTO sRunData ',...
     ' (singleid,masterid,cellid,batch,respath,resfile,kernfile)',...
     ' VALUES (',num2str(filedata.singleid),',',...
     num2str(filedata.masterid),',',...
     '"',filedata.cellid,'",',...
     num2str(batchid),',',...
     '"',respath,'",',...
     '"',[filedata.cellid,resbase],'",',...
     '"',[filedata.cellid,kernbase],'")'];
[aa,bb,runid]=mysql(sql);

fprintf('batchid=%d added runid %d for cell %s\n',...
        batchid,runid,cellid);


