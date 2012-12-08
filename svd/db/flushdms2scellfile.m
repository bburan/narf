% function [cellfileids]=checkraw2cell(rawids,fmt)
%
% (if not already done) convert to pype files to matlab-readable
% format and translate entries from gDataRaw to sCellFile
%
% if the conversion has already been done, simply return the
% cellfileids corresponding to the rawids and fmt (default='pixel')
% specified.
%
% CREATED SVD 3/3/03  (ripped off dbraw2scellfile.m)
%
function [cellfileids]=flushdms2scellfile(rawids,fmt)

dbopen;

% stimulus imsm archive path may vary with the well.
% after the end of a session, path should be updated to show
% archived locations of data.
stimpathdata.reese{13}='/auto/data/archive/stimarchive/REESE13/imsm/';
stimpathdata.reese{14}='/auto/data/archive/stimarchive/REESE14/imsm/';
stimpathdata.reese{15}='/auto/data/archive/stimarchive/REESE15/imsm/';
stimpathdata.reese{16}='/auto/data/archive/stimarchive/REESE16/imsm/';
stimpathdata.reese{17}='/auto/data/archive/stimarchive/REESE17/imsm/';
stimpathdata.mac{11}='/auto/data/archive/stimarchive/MAC11/imsm/';
stimpathdata.van{1}='/auto/data/archive/stimarchive/VAN01/imsm/';

% root of path to store processed data
% ie, data goes to:  $resproot/$pen_id
resproot='/auto/data/archive/reese_data/';

rcsetstrings;
if ~exist('fmt','var'),
   fmt='pixel';
end
stimfmtcode=find(strcmp(stimfilefmtstr,fmt))-1;
fprintf('stimfilefmt=%s (code=%d)\n',fmt,stimfmtcode);

RUNCLASSID=2;
rejectcells='("R240A","R247B","model","models")';

if ~exist('rawids'),
   GENERALUPDATE=1;
   sql=['SELECT * FROM gDataRaw',...
        ' WHERE runclassid=',num2str(RUNCLASSID),...
        ' AND bad=0',...
        ' AND not(cellid in ',rejectcells,')',...
        ' ORDER BY cellid,id'];
   rawdata=mysql(sql);
else
   GENERALUPDATE=0;
   srawid='(';
   for ii=1:length(rawids),
      srawid=[srawid,num2str(rawids(ii)),','];
   end
   srawid(end)=')';
   
   sql=['SELECT * FROM gDataRaw',...
        ' WHERE id in ',srawid,...
        ' AND bad=0',...
        ' ORDER BY cellid,id'];
   rawdata=mysql(sql);
end

rawids=cat(1,rawdata.id);
cellfileids=[];
%keyboard
for rawidx=1:length(rawdata),
   fprintf('** cellid %s rawid %d\n',...
           rawdata(rawidx).cellid,rawdata(rawidx).id);
   
   sql=['SELECT * FROM gCellMaster where id=',...
        num2str(rawdata(rawidx).masterid)];
   celldata=mysql(sql);
   
   sql=['SELECT * FROM sCellFile',...
        ' WHERE rawid=',num2str(rawdata(rawidx).id),...
        ' AND stimfilefmt="',fmt,'"'];
   cellfiledata=mysql(sql);
   
   if strcmp(celldata.animal,'reese') & celldata.well>13,
      fprintf('reese well>13: Redirecting stimpath & resppath!\n');
      
      stimpathidx=celldata.well;
      stimpathlist=getfield(stimpathdata,celldata.animal);
      
      if stimpathidx<=length(stimpathlist),
         stimpath=stimpathlist{stimpathidx};
         stimpath=stimpath(1:end-5);
      else
         fprintf('stimpath not found for %s!\n',stimfile{1});
         input('enter path: ',stimfile);
      end
      
      tsp=rawdata(rawidx).stimpath;
      trp=rawdata(rawidx).resppath;
      if strcmp(tsp(1:11),'/auto/stim/'),
         rawdata(rawidx).stimpath=[stimpath,tsp(12:end)];
      end
      if strcmp(trp(1:17),'/auto/reese_data/'),
         rawdata(rawidx).resppath=...
             ['/auto/data/archive/reese_data/',trp(18:end)];
      end
   end
   
   if length(cellfiledata)>0,
      fprintf('Has sCellFile entry stim=%s. Skipping...\n',...
              cellfiledata.stimfile);
      cellfileids=[cellfileids;cellfiledata.id];
      
      if ~exist([cellfiledata.stimpath,cellfiledata.stimfile],'file'),
         disp('stimfile not found!');
         keyboard;
      end
      if ~exist([cellfiledata.path,cellfiledata.respfile],'file'),
         disp('stimfile not found!');
         keyboard;
      end
      
   elseif ~strcmp(celldata.info(1:4),'CELL'),
      fprintf('No sCellFile entry but SGI-aged. Skipping...\n');
      
   elseif RUNCLASSID==1 & (strcmp(celldata.animal,'mac') | ...
          (strcmp(celldata.animal,'reese') & rawdata(rawidx).id<758)),
      disp('skipping data because old gratrev');
      
   else
      % is there any data in db from this rawdata entry?
      sql=['SELECT * FROM sCellFile',...
           ' WHERE rawid=',num2str(rawdata(rawidx).id)];
      cellfiledata=mysql(sql);
      
      % construct path for processed pype data
      penid=dbget('gCellMaster','penid',rawdata(rawidx).masterid);
      penname=dbget('gPenetration','penname',penid);
      resppath=[resproot,penname,'/'];
      if ~exist(resppath,'dir'),
         unix(['mkdir ',resppath]);
      end
      if length(cellfiledata)==0,
         % no existing entries. extract!
         matfile=processpypedata([rawdata(rawidx).resppath,...
                    rawdata(rawidx).respfile]);
         [respfile,oldpath]=basename(matfile);
      else
         respfile=cellfiledata(1).respfile;
         oldpath=cellfiledata(1).path;
         matfile=[oldpath,respfile];
      end
      
      if ~strcmp(oldpath,resppath),
         fprintf('Copying respfile to %s\n',resppath);
         unix(['cp ',matfile,' ',resppath]);
      end
      
      % figure out construction of pre-computed imsm file:
      stimbase=rawdata(rawidx).stimpath;
      if stimbase(end)=='/',
         stimbase=stimbase(1:(end-1));
      end
      ii=max(findstr(stimbase,'/'));
      jj=max([0 findstr(stimbase(1:ii-1),'/')]);
      ostimbase=stimbase((jj+1):end);
      
      % for new natrev naming scheme where dir is only size
      if ~isnan(str2double(stimbase((ii+1):end)))
         stimbase(ii)='-';
         ii=max(findstr(stimbase,'/'));
      end
      stimbase=[stimbase((ii+1):end),'.',rawdata(rawidx).stimfile,'.'];
      
      stimfile=[stimbase,fmt,'.imsm'];
      
      stimpathidx=celldata.well;
      stimpathlist=getfield(stimpathdata,celldata.animal);
      
      if stimpathidx<=length(stimpathlist),
         stimpath=stimpathlist{stimpathidx};
      else
         fprintf('stimpath not found for %s!\n',stimfile{1});
         input('enter path: ',stimfile);
      end
      
      % figure out location of pgms
      
      if ~exist([stimpath,stimfile],'file'),
         ii=findstr(stimpath,'imsm/');
         indexpath=[stimpath(1:ii-1),ostimbase,'/' ];
         indexfn=rawdata(rawidx).stimfile;
         fprintf('creating full pix imsm: %s%s\n',...
                 stimpath,stimfile);
         [framecount,iconside]=pypestim2imsmraw(indexpath,indexfn,...
                                                [stimpath,stimfile],1);
      else
         [framecount,iconside]=imfileinfo([stimpath,stimfile]);
         
      end
      
      % this following thing SHOULD work, but above hack is passable
      stimwindowsize=iconside(1);
      stimiconside=sprintf('%d,%d',iconside);
      
      z=load(matfile);
      resplen=sum(~isnan(z.psth(:,1)));
      respvarname='psth';
      respfmtcode=0;
      respfilefmt='PSTH'; % rather than PFTH
      
      stimspeedid=rawdata(rawidx).stimspeedid;
      
      [aff,rawdata(rawidx).cellfileid]=...
          sqlinsert('sCellFile',...
                    'cellid',rawdata(rawidx).cellid,...
                    'masterid',rawdata(rawidx).masterid,...
                    'rawid',rawdata(rawidx).id,...
                    'runclassid',rawdata(rawidx).runclassid,...
                    'path',resppath,...
                    'resplen',resplen,...
                    'respfile',respfile,...
                    'respvarname',respvarname,...
                    'respfiletype',1,...
                    'respfilefmt',respfilefmt,...
                    'respfmtcode',respfmtcode,...
                    'stimfile',stimfile,...
                    'stimpath',stimpath,...
                    'stimwindowsize',stimwindowsize,...
                    'stimfilefmt',fmt,...
                    'stimfmtcode',stimfmtcode,...
                    'stimspeedid',stimspeedid,...
                    'stimiconside',stimiconside,...
                    'addedby','david',...
                    'info','checkraw2cell.m');
      cellfileids=[cellfileids;rawdata(rawidx).cellfileid];
      
      % update rawdata stuff if not already there
      if isempty(rawdata(rawidx).stimfile),
         stimpath=z.h.play_dir0;
         if ~strcmp(stimpath(end),'/'),
            stimpath=[stimpath,'/'];
         end
         stimfile=z.h.index0;
         sql=['UPDATE gDataRaw SET',...
              ' stimpath="',stimpath,'",',...
              ' stimfile="',stimfile,'"',...
              ' WHERE id=',num2str(rawdata(rawidx).id)];
         mysql(sql);
      end
   end
end

