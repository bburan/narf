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

[cellfiledata,cellids]=dbbatchcells(batchid);

rejectcells={'R240A','R247B','modelss'};
cellids=setdiff(cellids,rejectcells);


if exist('matchbatchid','var') && ~isempty(matchbatchid),
   sql=['SELECT cellid FROM sRunData where batch=',num2str(matchbatchid)];
   matchdata=mysql(sql);
   matchcellid={matchdata.cellid};

   cellids=intersect(cellids,matchcellid);
end

rcsetstrings;

% switch from jlg to nsl labs
if batchid<108,
   respath=['/auto/k5/david/data/batch',num2str(batchid),'/'];
else
   respath=['/auto/data/nsl/users/svd/data/batch',num2str(batchid),'/'];
end

% find current entries in sRunData that don't have any sCellFiles
sql=['SELECT * FROM sRunData WHERE batch=',num2str(batchid)];
rundata=mysql(sql);

for ii=1:length(rundata),
   scdata=dbbatchcells(batchid,rundata(ii).cellid);
   if length(scdata)==0,
      fprintf('cell %s missing... deleting sRunData entry\n',...
              rundata(ii).cellid);
      yn=input('Delete (y/[n])? ','s');
      if ~isempty(yn) && strcmp(yn(1),'y'),
      
         sql=['DELETE FROM sRunData WHERE id=', ...
              num2str(rundata(ii).id)];
         mysql(sql);
         sql=['DELETE FROM sResults WHERE runid=', ...
              num2str(rundata(ii).id)];
         mysql(sql);
      end
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
for ii=1:length(cellids),
   
   sql=['SELECT * FROM sRunData WHERE cellid="',cellids{ii},'"',...
        ' AND batch=',num2str(batchid)];
   rundata=mysql(sql);
   
   if length(rundata)>0,
      %fprintf('Entry exists for cellid %s. Skipping.\n',cellids{ii});
      
   elseif (strcmp(cellids{ii}(1:5),'sir00') || ...
         strcmp(cellids{ii}(1:5),'sir01') || ...
         strcmp(cellids{ii}(1:5),'sir02') || ...
         strcmp(cellids{ii}(1:5),'sir03') || ...
         strcmp(cellids{ii}(1:5),'sir04') || ...
         strcmp(cellids{ii}(1:5),'sir05')) && ...
         ismember(batchid,[133:138 185:188 194]),
      
      fprintf('Skipping cellid %s for technical reason.\n',cellids{ii});
         
   else
      if yestoall,
         fprintf('Auto-adding entry for cellid %s\n',cellids{ii});
         yn='y';
      else
         yn=input(sprintf('Add entry for cellid %s (y/[n]) ',...
                          cellids{ii}),'s');
      end
      if ~isempty(yn) && strcmp(yn(1),'y'),
         sql=['SELECT * FROM gSingleCell WHERE cellid="',cellids{ii},'"'];
         singledata=mysql(sql);
         if isempty(singledata),
            clear singledata
            singledata.id=0;
            singledata.masterid=0;
            warning('gSingleCell entry missing.');
         end
         %keyboard
         sql=['INSERT INTO sRunData ',...
              ' (singleid,masterid,cellid,batch,respath,resfile,kernfile)',...
              ' VALUES (',num2str(singledata.id),',',...
              num2str(singledata.masterid),',',...
              '"',cellids{ii},'",',...
              num2str(batchid),',',...
              '"',respath,'",',...
              '"',[cellids{ii},resbase],'",',...
              '"',[cellids{ii},kernbase],'")'];
         [aa,bb,newrunid]=mysql(sql);
         runids=[runids;newrunid];
         runcount=runcount+1;
      end
   end
end

fprintf('batchid=%d added %d runs, resbase=%s\n',...
        batchid,runcount,resbase);
