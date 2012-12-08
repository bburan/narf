
dbopen;
sql=['SELECT gSingleCell.*,gPenetration.numchans FROM gSingleCell',...
     ' INNER JOIN gPenetration ON gSingleCell.penid=' ...
     ' gPenetration.id WHERE not(gPenetration.training)'];
singledata=mysql(sql);

for ii=1:length(singledata),
   sid=singledata(ii).id;
   ocellid=singledata(ii).cellid;
   siteid=singledata(ii).siteid;
   channum=singledata(ii).channum;
   unit=singledata(ii).unit;
   if unit==0,
      unit=str2num(ocellid(end));
      if isempty(unit),
         unit=1;
      end
   end
   if singledata(ii).numchans<=8,
      cellid=[siteid '-' char('a'-1+channum) num2str(unit)];
   else
      cellid=[siteid '-' sprintf('%02d',channum) '-' ...
              num2str(unit)];
   end
   if ~strcmp(cellid,ocellid),
      fprintf('%d/%d. cellid mismatch: %s v %s\n',ii,sid,cellid,ocellid);
      [res,aff]=mysql(['UPDATE gSingleCell set cellid="',cellid,'"',...
                       ' WHERE id=',num2str(sid)]);
      fprintf('gSingleCell: %d rows updated\n',aff);
      [res,aff]=mysql(['UPDATE gSingleRaw set cellid="',cellid,'"',...
           ' WHERE cellid="',ocellid,'"']);
      fprintf('gSingleRaw: %d rows updated\n',aff);
      [res,aff]=mysql(['UPDATE sCellFile set cellid="',cellid,'"',...
           ' WHERE cellid="',ocellid,'"']);
      fprintf('sCellFile: %d rows updated\n',aff);
      [res,aff]=mysql(['UPDATE sRunData set cellid="',cellid,'"',...
           ' WHERE cellid="',ocellid,'"']);
      fprintf('sRunData: %d rows updated\n',aff);
   end
   if unit~=singledata(ii).unit,
      fprintf('%d/%d. %s unit mismatch: %d v %d\n',ii,sid,cellid,unit,...
              singledata(ii).unit);
      [res,aff]=mysql(['UPDATE gSingleCell set unit=',num2str(unit),...
                       ' WHERE id=',num2str(sid)]);
      fprintf('gSingleCell: %d rows updated\n',aff);
      [res,aff]=mysql(['UPDATE gSingleRaw set unit=',num2str(unit),...
                       ' WHERE singleid=',num2str(sid)]);
      fprintf('gSingleRaw: %d rows updated\n',aff);
      
   end
end


cd /auto/p1/svd/docs/current/lfp_vs_stim/

load TORData

dbopen;
b=178;
for ii=1:length(TORfname),
   sql=['SELECT min(cellid) as cellid,channum,masterid',...
        ' FROM sCellFile where respfile="',TORfname{ii},'"',...
        ' GROUP BY channum'];
   cfdata=mysql(sql);
   
   for cc=1:length(cfdata),
      sql=['SELECT * FROM sRunData WHERE batch=',num2str(b),...
           ' AND cellid="',cfdata(cc).cellid,'"'];
      rdata=mysql(sql);
      if isempty(rdata),
         respath=sprintf('/auto/data/nsl/users/svd/data/batch%d/',b);
         resfile=sprintf('%s.TOR.LFP.%d.res.mat',...
                         cfdata(cc).cellid,b);
         
         sqlinsert('sRunData',...
                     'masterid',cfdata(cc).masterid,...
                     'cellid',cfdata(cc).cellid,...
                     'batch',b,...
                     'respath',respath,...
                     'resfile',resfile);
         fprintf('added %s: %s\n',cfdata(cc).cellid,resfile);
      end        
   end
end



baphy_set_path;


dbopen;

sql=['SELECT * FROM gData WHERE name="Ref_NoiseType"',...
     ' AND svalue="SpectSmooth"'];
datadata=mysql(sql);

for ii=1:length(datadata),
   sql=['UPDATE gData SET value=0 WHERE name="Ref_SNR"',...
        ' AND rawid=',num2str(datadata(ii).rawid)];
   mysql(sql);
   
   sql=['UPDATE sCellFile SET stimsnr=0',...
        ' WHERE rawid=',num2str(datadata(ii).rawid)];
   mysql(sql);
   
end


sql='select * from sRunData where id>32604;';
rdata=mysql(sql);
for ii=1:length(rdata),
   isrmakequeue(rdata(ii).cellid,rdata(ii).batch,'cellxcmaster');
end




sql=['SELECT gDataRaw.* FROM gDataRaw LEFT JOIN gData',...
     ' ON (gData.rawid=gDataRaw.id AND gData.name="DI_err")',...
     ' WHERE runclass="DMS" and not(bad)',...
     ' AND gDataRaw.id>=28218',...
     ' AND isnull(gData.value)'];

rawfiledata=mysql(sql);

for ii=1:length(rawfiledata),
   mfile=[rawfiledata(ii).resppath rawfiledata(ii).parmfile];
   fprintf('\n%d/%d: %s\n',ii,length(rawfiledata),basename(mfile));
   if exist(mfile,'file'),
      dms_plot_print(mfile);
      close(gcf);
   end
end


return

d=load('/auto/p1/svd/code/autils/TNTDataforStephen.mat');
sitelist={};
for ii=1:length(d.name),
   if ~sum(strcmp(sitelist,d.name{ii}(1:5))),
      sitelist{length(sitelist)+1}=d.name{ii}(1:5);
   end
end

dbopen
sql=['SELECT * FROM gCellMaster WHERE siteid like "v0%" AND siteid<="v042"' ...
     ' AND (isnull(area) or area="") AND not(training)'];
sitedata=mysql(sql);

for ii=1:length(sitedata),
   
      sql=['SELECT * FROM gPenetration WHERE id = ',num2str(sitedata(ii).penid)];
      pendata=mysql(sql);
      
      masterid=sitedata(ii).id;
      
      area='A1';
      areastr='';
      for jj=1:pendata.numchans,
         areastr=[areastr,area,','];
      end
      areastr=areastr(1:(end-1));
      
      sql=['UPDATE gCellMaster set area="',areastr,'"',...
           ' WHERE penid=',num2str(pendata.id)];
      mysql(sql);
      
      sql=['UPDATE gSingleCell set area="',area,'"',...
           ' WHERE penid=',num2str(pendata.id)];
      mysql(sql);
      
      
      sql=['UPDATE sCellFile set area="',area,'"',...
           ' WHERE masterid=',num2str(masterid)];
      mysql(sql);
end

return

dbopen
sql=['SELECT * FROM gCellMaster WHERE siteid like "mer%" AND not(isnull(area)) AND area<>""'];
sitedata=mysql(sql);

for ii=1:length(sitedata),
   revisedarea=strrep(sitedata(ii).area,'pac','a1');
   revisedarea=strrep(revisedarea,'pfc','fc');
   if ~strcmp(revisedarea,sitedata(ii).area),
      revisedarea
      sql=['UPDATE gCellMaster set area="',revisedarea,'"',...
           ' WHERE id=',num2str(sitedata(ii).id)];
      mysql(sql);
   end
   
   areas=strsep(revisedarea,',');
   for jj=1:length(areas),
      area=areas{jj};
      
      sql=['UPDATE gSingleCell set area="',area,'"',...
           ' WHERE masterid=',num2str(sitedata(ii).id),...
              ' AND channum=',num2str(jj)];
      mysql(sql);
      
      sql=['UPDATE sCellFile set area="',area,'"',...
           ' WHERE masterid=',num2str(sitedata(ii).id),...
           ' AND channum=',num2str(jj)];
      mysql(sql);
   end
end

return

dbopen
sql=['SELECT cellid,count(id) as filecount FROM sCellFile',...
     ' WHERE (cellid LIKE  "plu2%"',...
     ' OR cellid LIKE  "plu3%"',...
     ' OR cellid LIKE  "a3%")',...
     ' AND respfile like "%_a_%"',...
     ' GROUP BY cellid ORDER BY cellid'];
celldata=mysql(sql);

badlist=[];
for ii=130:length(celldata),
   try
      ii
      close all
      quicksum(celldata(ii).cellid);
      drawnow;
      print -f1 -Poud
   catch
      badlist=[badlist ii]
   end
   
end
return


dbopen

sql=['SELECT * FROM gPenetration WHERE not(training)',...
     ' AND animal="pluto" AND penname like "plu4%"'];
pendata=mysql(sql);

for ii=1:length(pendata),
   area='FC';
   areastr='';
   for jj=1:pendata(ii).numchans,
      areastr=[areastr,area,','];
   end
   areastr=areastr(1:(end-1));
   
   sql=['UPDATE gCellMaster set area="',areastr,'"',...
        ' WHERE penid=',num2str(pendata(ii).id)];
   mysql(sql);
   
   sql=['UPDATE gSingleCell set area="',area,'"',...
        ' WHERE penid=',num2str(pendata(ii).id)];
   mysql(sql);
   
   sql=['SELECT distinct masterid FROM gSingleCell',...
        ' WHERE  penid=',num2str(pendata(ii).id)];
   singledata=mysql(sql);
   
   for jj=1:length(singledata),
      masterid=singledata(jj).masterid;
      
      sql=['UPDATE sCellFile set area="',area,'"',...
           ' WHERE masterid=',num2str(masterid)];
      mysql(sql);
   end
end

return

   
sql=['SELECT cellid,count(id) as filecount FROM sCellFile',...
     ' WHERE (cellid LIKE  "plu2%"',...
     ' OR cellid LIKE  "plu3%"',...
     ' OR cellid LIKE  "a3%")',...
     ' AND respfile like "%_a_%"',...
     ' GROUP BY cellid ORDER BY cellid'];
celldata=mysql(sql);

badlist=[];
for ii=130:length(celldata),
   try
      ii
      close all
      quicksum(celldata(ii).cellid);
      drawnow;
      print -f1 -Poud
   catch
      badlist=[badlist ii]
   end
   
end
return


dbopen;
sql=['SELECT gDataRaw.* FROM gDataRaw LEFT JOIN gData',...
     ' ON (gData.rawid=gDataRaw.id AND gData.name="Overlay_Ref_Tar")',...
     ' WHERE runclass="DMS" and not(bad) and not(training)',...
    ' AND isnull(gData.value)'];

rawfiledata=mysql(sql);

for ii=1:length(rawfiledata),
   ii
   if exist([rawfiledata(ii).resppath rawfiledata(ii).parmfile],'file'),
      clear exptparams globalparams exptevents events
      LoadMFile([rawfiledata(ii).resppath rawfiledata(ii).parmfile]);
      
      Parameters=[];
      if isfield(exptparams,'targ_atten'),
         Parameters.Tar_Atten=exptparams.targ_atten;
      end
      if isfield(exptparams,'targ_atten'),
         Parameters.Overlay_Ref_Tar=exptparams.overlay_reftar;
      end
      if isfield(exptparams,'targ_atten'),
         Parameters.Use_Catch=exptparams.use_catch;
      end
      if ~isempty(Parameters),
         rawfiledata(ii).parmfile
         Parameters
         dbWriteData(globalparams.rawid,Parameters,0,1);
      end
   end
end 

return

dbopen;

sql=['SELECT gSingleRaw.*,gDataRaw.matlabfile FROM gSingleRaw,gDataRaw',...
     ' WHERE gSingleRaw.rawid=gDataRaw.id',...
     ' AND (gSingleRaw.isolation=0 OR isnull(gSingleRaw.isolation))',...
     ' AND gSingleRaw.info<>"populatensl.m"',...
     ' AND gDataRaw.matlabfile<>""',...
     ' AND not(isnull(gDataRaw.matlabfile))'];

rawfiledata=mysql(sql);

for ii=189:length(rawfiledata),
   
   channum=rawfiledata(ii).channum;
   unit=rawfiledata(ii).unit;
   
   p=load(rawfiledata(ii).matlabfile);
   
   if isfield(p,'sortextras') && ...
         length(p.sortextras)>=channum && ...
         ~isempty(p.sortextras{channum}) && ...
         p.sortinfo{channum}{1}(1).Ncl>=unit,
      if isfield(p.sortextras{channum},'unitmean'),
         snr=std(p.sortextras{channum}.unitmean(:,unit))./...
             p.sortextras{channum}.sigma;
      else
         snr=std(mean(p.sortinfo{channum}{1}(unit).env(2:3,:),1))./...
             p.sortextras{channum}.sigma;
         
      end
      if ~isnan(snr),
         isopct=round( 100*erf(snr./2) .* 10) ./10
         if isopct==0,
            isopct=-1;  % tag as not valid
         end
         
         fprintf('%d: %s, %s: %.1f\n',ii,basename(rawfiledata(ii).matlabfile),...
                 rawfiledata(ii).cellid,isopct);
         
         sql=['UPDATE gSingleRaw set isolation=',num2str(isopct),...
              ' WHERE id=',num2str(rawfiledata(ii).id)];
         mysql(sql);
      end
   end
   
end
   

return



dbopen;

sql=['SELECT gData.value,sCellFile.* FROM sCellFile LEFT JOIN gData',...
     ' ON sCellFile.rawid=gData.rawid AND gData.name="Ref_SNR"'];
cellfiledata=mysql(sql);

for ii=1:length(cellfiledata),
   fprintf('%d %s\n',ii,cellfiledata(ii).respfile);
   if isempty(cellfiledata(ii).value),
      if cellfiledata(ii).stimsnr~=1000,
         sql=['UPDATE sCellFile set stimsnr=1000 WHERE id=',num2str(cellfiledata(ii).id)]
         mysql(sql);
      end
   else
      sql=['UPDATE sCellFile set stimsnr=',num2str(cellfiledata(ii).value),...
           ' WHERE id=',num2str(cellfiledata(ii).id)]
      mysql(sql);
   end

end

return

dbopen;

sql=['SELECT gData.value,sCellFile.* FROM sCellFile LEFT JOIN gData',...
     ' ON sCellFile.rawid=gData.rawid AND gData.name="Ref_Duration"'];
cellfiledata=mysql(sql);

for ii=1:length(cellfiledata),
   fprintf('%d %s\n',ii,cellfiledata(ii).respfile);
   if isempty(cellfiledata(ii).value),
      sql=['UPDATE sCellFile set stimspeedid=0 WHERE id=',num2str(cellfiledata(ii).id)]
   else
      sql=['UPDATE sCellFile set stimspeedid=',num2str(cellfiledata(ii).value),...
           ' WHERE id=',num2str(cellfiledata(ii).id)]
   end 
   mysql(sql);
end


dbopen;

sql=['SELECT * FROM gDataRaw where not(training) AND not(bad) and not(cellid like "tst%") and not(cellid like "test%")'];
rawdata=mysql(sql);

parmmissinglist=[];
rawmissinglist=[];
evpmissinglist=[];
for ii=1:length(rawdata),
   [parmfile,resppath]=basename(rawdata(ii).parmfile);
   
   if isempty(resppath),
      resppath=rawdata(ii).resppath;
   end
   
   if resppath(2)~=':', % is, not a dos path
      [respfileevp,tt]=basename(rawdata(ii).respfileevp);
      if ~isempty(tt) & ~strcmp(tt,resppath),
         disp('parmpath evppath mismatch!');
         keyboard
      end
      if isempty(respfileevp),
         [p,b,e]=fileparts(parmfile);
         respfileevp=[b '.evp'];
      end
      matlabfile=rawdata(ii).matlabfile;
      
      
      if (~isempty(matlabfile) & ~exist(matlabfile,'file')),
         [resppath matlabfile]
         disp('spike file missing');
      end
      if ~exist([resppath parmfile],'file'),
         [resppath parmfile]
         disp('parm file missing');
         parmmissinglist=[parmmissinglist ii];
      end
      
      if (~exist([resppath respfileevp],'file') & ...
          ~exist([resppath respfileevp '.gz'],'file')),
         [resppath respfileevp]
         disp('evp file missing');
         evpmissinglist=[evpmissinglist ii];
      end
   end
end

return


dbopen;

sql=['SELECT * FROM gDataRaw where not(bad) AND parmfile like "/afs/%"'];
rawdata=mysql(sql);

oldroot='/afs/glue.umd.edu/department/isr/labs/nsl/projects/daqsc/';
newroot='/auto/data/daq/';
oldroot2='/afs/glue.umd.edu/department/isr/labs/nsl/';
newroot2='/auto/data/nsl/';
rawmissinglist=[];
evpmissinglist=[];
for ii=1:length(rawdata),
   ii
   [parmfile,resppath]=basename(rawdata(ii).parmfile);
   if ~isempty(rawdata(ii).resppath),
      rawdatapath=rawdata(ii).resppath;
      respfile=[rawdatapath rawdata(ii).respfile];
   else
      respfile=rawdata(ii).respfile;
   end
   
   if isempty(rawdata(ii).respfileevp),
      [p,b,e]=fileparts(parmfile);
      respfileevp=[b '.evp'];
      
   elseif strcmp(rawdata(ii).respfileevp(1:5),'/afs/'),
      [respfileevp,tt]=basename(rawdata(ii).respfileevp);
      if ~strcmp(tt,resppath),
         disp('parmpath evppath mismatch!');
         keyboard
      end
   end
   matlabfile=rawdata(ii).matlabfile;
   
   resppath=strrep(resppath,oldroot,newroot);
   respfile=strrep(respfile,oldroot,newroot);
   matlabfile=strrep(matlabfile,oldroot,newroot);
   matlabfile=strrep(matlabfile,oldroot2,newroot2);
   [resppath parmfile]
   
   if ~exist([resppath parmfile],'file') | ...
         (~isempty(matlabfile) & ~exist(matlabfile,'file')),
      disp('file missing');
      keyboard
   end
   if (~exist([resppath respfileevp],'file') & ...
          ~exist([resppath respfileevp '.gz'],'file')),
      disp('evp file missing');
      evpmissinglist=[evpmissinglist ii];
   end
   if ~exist(respfile,'file'),
      disp('raw file missing!');
      rawmissinglist=[rawmissinglist ii];
   end
   
   sql=['UPDATE gDataRaw set resppath="',resppath,'",'...
        ' respfile="',respfile,'",'...
        ' matlabfile="',matlabfile,'",'...
        ' parmfile="',parmfile,'",',...
        ' respfileevp="',respfileevp,'"'...
        ' WHERE id=',num2str(rawdata(ii).id)]
   mysql(sql);
end

return


dbopen;

sql=['SELECT sCellFile.*,gDataRaw.resppath as rawpath,parmfile',...
     ' FROM sCellFile,gDataRaw WHERE sCellFile.rawid=gDataRaw.id',...
     ' AND sCellFile.stimpath like "/auto/data/%"'];
cellfiledata=mysql(sql);

for ii=1:length(cellfiledata),
   if cellfiledata(ii).parmfile(1)=='/' | ...
          cellfiledata(ii).parmfile(2)==':',
      parmfile=cellfiledata(ii).parmfile;
      [bb,pp]=basename(cellfiledata(ii).parmfile);
   else
      bb=cellfiledata(ii).parmfile;
      pp=cellfiledata(ii).rawpath;
   end
   pp=strrep(pp,'/afs/glue.umd.edu/department/isr/labs/nsl/projects/daqsc',...
             '/auto/data/daq');
   pp=strrep(pp,'M:/daq',...
             '/auto/data/daq');
   if ~strcmp(bb,cellfiledata(ii).stimfile) & ...
         strcmp(bb,[cellfiledata(ii).stimfile '.m']),
      cellfiledata(ii).stimfile=bb;
   end
   if ~strcmp(bb,cellfiledata(ii).stimfile),
      disp('mismatch!!!');
      ii
   end
   
   if ~strcmp(pp,cellfiledata(ii).stimpath),
      
      fprintf('%d %s %s vs %s %s\n',cellfiledata(ii).rawid,...
              cellfiledata(ii).stimpath,cellfiledata(ii).stimfile,...
              pp,bb);
      sql=['UPDATE sCellFile set stimpath="',pp,'"',...
           ' WHERE id=',num2str(cellfiledata(ii).id)]
      mysql(sql);
   end
end
return


dbopen

sql=['SELECT * FROM gPenetration'];
pendata=mysql(sql);

for ii=1:length(pendata),
   if pendata(ii).pendate(3)=='/',
      t=strsep(pendata(ii).pendate,'/',1);
      newdate=[t{3} '-' t{1} '-' t{2}];
      fprintf('%s -> %s\n',pendata(ii).pendate,newdate);
   elseif pendata(ii).pendate(3)=='-',
      t=datenum(pendata(ii).pendate);
      newdate=datestr(t,'yyyy-mm-dd');
      fprintf('%s -> %s\n',pendata(ii).pendate,newdate);
   elseif pendata(ii).pendate(5)=='-',
      fprintf('already ok: %s\n',pendata(ii).pendate);
   end
   sql=['UPDATE gPenetration set pendate="' newdate ...
        '" WHERE id=' num2str(pendata(ii).id)];
   mysql(sql);
end


return


dbopen

sql=['SELECT DISTINCT sRunData.id,sCellFile.singleid',...
     ' FROM sRunData INNER JOIN sCellFile',...
     ' ON sRunData.cellid=sCellFile.cellid',...
     ' WHERE sRunData.singleid=0'];
rundata=mysql(sql);
sql=['SELECT DISTINCT sRunData.id,sCellFile.singleid',...
     ' FROM sRunData INNER JOIN sCellFile',...
     ' ON sRunData.cellid=sCellFile.cellid',...
     ' WHERE sRunData.singleid=0'];
rundata=mysql(sql);

for ii=1:length(rundata),
   sql=['UPDATE sRunData set singleid=',num2str(rundata(ii).singleid),...
        ' WHERE id=',num2str(rundata(ii).id)];
   mysql(sql);
end



dbopen

for batchid=87:107,
   sql=['SELECT * FROM sRunData WHERE batch=',num2str(batchid)];
   rundata=mysql(sql);
   if length(rundata)>0,
      resfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
      
      if ~exist(resfile,'file'),
         resfile=[rundata(1).respath,rundata(1).kernfile,'.gz'];
      end
      
      if exist(resfile,'file'),
         z=zload(resfile);
         
         if isfield(z,'params'),
            if isfield(z.params,'parmstring'),
               parmstring=char(z.params.parmstring);
            else
               parmstring='';
            end
            sql=['UPDATE sBatch set parmstring="',parmstring,...
                 '" WHERE id=',num2str(batchid)]
            mysql(sql);
            
            if 0
            stimfiltstr='';
            for ii=1:length(z.params.stimfilterparms),
               if isnumeric(z.params.stimfilterparms{ii}),
                  stimfiltstr=[stimfiltstr,',',...
                               num2str(z.params.stimfilterparms{ii})];
               else
                  stimfiltstr=[stimfiltstr,',''',...
                               z.params.stimfilterparms{ii},''''];
               end
            end
            if length(z.params.stimfilterparms)>0,
               stimfiltstr=stimfiltstr(2:end);
            end
            sql=['UPDATE sBatch set stimfilterparms="',stimfiltstr,...
                 '" WHERE id=',num2str(batchid)]
            mysql(sql);
            end
         end
      end
   end
end

return
   

dbopen;

sql=['SELECT * FROM sCellFile WHERE singleid=0 OR singlerawid=0'];
filedata=mysql(sql);

for ii=1:length(filedata),
   if filedata(ii).singleid==0 | filedata(ii).singlerawid==0,
      sql=['SELECT * FROM gSingleRaw',...
           ' WHERE rawid=',num2str(filedata(ii).rawid),...
           ' AND cellid="',filedata(ii).cellid,'"'];
      singledata=mysql(sql);
      if length(singledata)==1,
         sql=['UPDATE sCellFile SET',...
              ' singleid=',num2str(singledata.singleid),',',...
              ' singlerawid=',num2str(singledata.id),...
              ' WHERE id=',num2str(filedata(ii).id)];
         mysql(sql);
      else
         fprintf('no singledata for %s rawid=%d\n',...
                 filedata(ii).cellid,filedata(ii).rawid);
      end
   end
end 


return


dbopen;

sql=['SELECT * FROM gCellMaster'];
masterdata=mysql(sql);

for ii=1:length(masterdata),
   sql=['SELECT * FROM gSingleCell WHERE masterid=',...
        num2str(masterdata(ii).id)];
   singledata=mysql(sql);
   if length(singledata)==0,
      
      hp=masterdata(ii).handplot;
      ff=find(hp=='"');
      for ffidx=length(ff):-1:1,
         hp=[hp(1:ff(ffidx)) hp(ff(ffidx):end)];
      end
      if isempty(masterdata(ii).xoffset),
         xoffset=0;
      else
         xoffset=masterdata(ii).xoffset;
      end
      if isempty(masterdata(ii).yoffset),
         yoffset=0;
      else
         yoffset=masterdata(ii).yoffset;
      end
      
      % assume siteid=cellid (no mult cells per site yet!)
      sql=['INSERT INTO gSingleCell (siteid,cellid,masterid,penid,channel,'...
           'unit,area,handplot,rfsource,rfsize,xoffset,yoffset,'...
           'quality,latency,addedby,info,crap) VALUES (',...
           '"',masterdata(ii).cellid,'",',...
           '"',masterdata(ii).cellid,'",',...
           num2str(masterdata(ii).id),',',...
           num2str(masterdata(ii).penid),',',...
           '"a",1,',...
           '"',masterdata(ii).area,'",',...
           '"',hp,'",',...
           '"',masterdata(ii).rfsource,'",',...
           num2str(max([masterdata(ii).rfsize 0])),',',...
           num2str(xoffset),',',...
           num2str(yoffset),',',...
           num2str(masterdata(ii).quality),',',...
           num2str(max([masterdata(ii).latency 0])),',',...
           '"',masterdata(ii).addedby,'",',...
           '"',masterdata(ii).info,'",',...
           num2str(masterdata(ii).crap),')',...
          ];
      [res,aff,singleid]=mysql(char(sql));
   else
      singleid=singledata(1).id;
   end
   
   sql=['SELECT * FROM gDataRaw WHERE masterid=',num2str(masterdata(ii).id)];
   rawdata=mysql(sql);
   
   for jj=1:length(rawdata),
      sql=['SELECT * FROM gSingleRaw',...
           ' WHERE singleid=',num2str(singleid),...
           ' AND rawid=',num2str(rawdata(jj).id)];
      singledata=mysql(sql);
      
      if length(singledata)==0
         % singleraw entry doesn't exist yet. create one.
         sqlinsert('gSingleRaw',...
                   'cellid',masterdata(ii).cellid,...
                   'masterid',masterdata(ii).id,...
                   'singleid',singleid,...
                   'penid',masterdata(ii).penid,...
                   'rawid',rawdata(jj).id,...
                   'channel','p',...
                   'unit','1',...
                   'crap',masterdata(ii).crap,...
                   'isolation',max([rawdata(jj).isolation -101]),...
                   'addedby',masterdata(ii).addedby,...
                   'info',masterdata(ii).info);
      end
   end
end

return



dbopen;

sql=['SELECT * FROM gDataRaw WHERE resppath like "%R_DATA_ARCHIVE%"' ...
     ' and respfile=""'];
rawdata=mysql(sql);

for ii=1:length(rawdata),
   
   jj=10;
   matchfile=[2 2];
   
   while length(matchfile)>1
      matchfile=jls([rawdata(ii).resppath,...
                     rawdata(ii).matlabfile(27:(27+jj)),'*.d']);
      jj=jj+1;
   end
   
   if length(matchfile)==0,
      fprintf('%d %s   %d matches\n matfile: %s\n respfile: %s\n',...
              ii,rawdata(ii).cellid,length(matchfile),...
              rawdata(ii).matlabfile(27:end),...
              'NO MATCH');
      
   else
      fprintf('%d %s   %d matches\n matfile: %s\n respfile: %s\n',...
              ii,rawdata(ii).cellid,length(matchfile),...
              rawdata(ii).matlabfile(27:end),...
              matchfile{1});
      respfile=matchfile{1}(52:end);
      
      q=sprintf('want to make respfile %s (y/[n])? ',respfile);
      yn=input(q,'s');
      if strcmp(yn,'y'),
         % change entry in gDataRaw
         sql=['UPDATE gDataRaw SET respfile="',respfile,'"',...
              ' WHERE id=',num2str(rawdata(ii).id)]
         %keyboard
         mysql(sql);
      end
   end
end


